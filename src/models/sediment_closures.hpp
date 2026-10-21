#include <cmath>
#include "hydraulics.hpp"

namespace sediment {

    inline constexpr double g = 9.81;
    inline constexpr double RHO_W = 1000.0;  // water density, kg/m^3
    inline constexpr double NU = 1e-6;  // kinematic viscosity, m^2/s
    inline constexpr double KAPPA = 0.4;  // von Karman constant

    struct BedSourcedFlux {
        double flux;  // suspended bed material transport rate, kg/s
        double rouse; // Z = w_s / kappa u_star (dimensionless)
        double c_ref;  // reference concentration
        bool suspending;  // false if the bed stress is below the critical level
    };

    inline double settlingVelocity(double d50_mm, double rho_s) {
        // Implementing standard Stokes settling velocity
        double rho_effective = rho_s - RHO_W;
        double mu = NU * RHO_W;
        double d50_m = d50_mm / 1000.0;
        double stokes_velocity = rho_effective * g * std::pow(d50_m, 2) / (18 * mu);
        return stokes_velocity;
    }

    inline double dStar(double d50_mm, double rho_s) {
        // Nondimensional grain size
        double s = rho_s / RHO_W;
        double d50_m = d50_mm / 1000.0;
        double d_star = std::pow(g * (s - 1) / std::pow(NU, 2), 1.0 / 3) * d50_m;
        return d_star;
    }

    inline double criticalShields(double d50_mm, double rho_s) {

        // // Implementing Soulsby & Whitehouse (1997) critical shields parameter expression
        // // This should only be applied to non-cohesive sediments, so be wary of results for
        // // d50 < 62.5 microns
        // double D_star = d_star(d50_mm, rho_s);
        //
        // // Equation 10
        // double shields_crit = 0.30 / (1 + 1.2 * D_star) + 0.055 * (1 - std::exp(-0.02 * D_star));
        // return shields_crit;

        // Implementing van Rijn 1984 for consistency with the pickup function
        double D_star = dStar(d50_mm, rho_s);
        double shields_crit;  // van Rijn (1984a) Shields fit
        if (D_star <= 4.0) {
            shields_crit = 0.24  / D_star;
        }
        else if (D_star <= 10.0) {
            shields_crit = 0.14  * std::pow(D_star, -0.64);
        }
        else if (D_star <= 20.0) {
            shields_crit = 0.04  * std::pow(D_star, -0.10);
        }
        else if (D_star <= 150.0) {
            shields_crit = 0.013 * std::pow(D_star,  0.29);
        }
        else {
            shields_crit = 0.055;
        }
        return shields_crit;
    }

    inline double vanRijnCref(double d50_mm, double rho_s, double h, double shields, double shields_crit) {
        double D_star = dStar(d50_mm, rho_s);
        double T = (shields - shields_crit) / shields_crit;
        double d50_m = d50_mm / 1000.0;
        double a = std::max(2.5 * d50_m, 0.01 * h);  // Nikuradse roughness ks or 1% depth
        double c_ref = 0.015 * (d50_m / a) * (std::pow(T, 1.5) / std::pow(D_star, 0.3));
        return c_ref;
    }

    inline BedSourcedFlux bedSourcedSuspension(const hydraulics::HydraulicState& hs, double B, double d50_mm, double rho_s) {
        // Suspended sediment flux

        BedSourcedFlux flux;

        double shields_crit = criticalShields(d50_mm, rho_s);
        flux.suspending = hs.shields > shields_crit;
        if (!flux.suspending) {
            return flux;
        }

        double ws = settlingVelocity(d50_mm, rho_s);


        // Parameters for log law velocity profile
        double d50_m = d50_mm / 1000.0;
        double z0 = 2.5 * d50_m / 30;

        // And for rouse profile
        double Z = ws / (KAPPA * hs.ustar);
        double c_ref = vanRijnCref(d50_mm, rho_s, hs.depth, hs.shields, shields_crit);
        double z_ref = std::max(2.5 * d50_m, 0.01 * hs.depth);

        int N = 100; // Vertical discretization
        double dz = (hs.depth - z_ref) / N;
        double c_flux = 0.0;
        for (int i = 0; i <= N; ++i) {
            double z = z_ref + i * dz;

            // velocity from log law
            double u = (hs.ustar / KAPPA) * std::log(z / z0);

            // Concentration from Rouse profile
            double c = rho_s * c_ref * std::pow((std::max(hs.depth - z, 0.0) / z) * (z_ref / (hs.depth - z_ref)), Z);

            double w = (i == 0 || i == N) ? 0.5 : 1.0;   // trapezoid weights
            c_flux += w * c * u * dz; // kg / m / s
        }

        flux.flux = c_flux * B; // multiply by width to get kg / s
        flux.rouse = Z;
        flux.c_ref = c_ref;

        return flux;
    }
} // namespace sediment