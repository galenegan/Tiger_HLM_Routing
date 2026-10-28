#pragma once
#include <algorithm>
#include <cmath>
#include "hydraulics.hpp"

// Suspended bed-material load following van Rijn (1984a, Part I: J. Hydraul. Eng. 110(10);
// 1984b, Part II: 110(11)).
// Deviations: Ferguson-Church fall velocity of D50 in place of Part II eqs. 11-13 and 39, plane-bed
// roughness 3 D90 for the reference level, and the log-law profile scaled to the Manning
// mean velocity. The grain-size functions assume d50_mm > 0 so bedSourcedSuspension enforces it.

namespace sediment {

    inline constexpr double g = 9.81;
    inline constexpr double RHO_W = 1000.0;  // water density, kg/m^3
    inline constexpr double NU = 1e-6;  // kinematic viscosity, m^2/s
    inline constexpr double KAPPA = 0.4;  // von Karman constant
    inline constexpr double USTAR_MIN = 1e-6;  // floor on u_star in the Rouse number, m/s

    // Ferguson and Church (2004) drag constants, natural sand grains
    inline constexpr double FC_C1 = 18.0;
    inline constexpr double FC_C2 = 1.0;

    // The per-link data carry D50 only
    inline constexpr double D90_OVER_D50 = 1.5;
    // Maximum (bed) volumetric concentration, van Rijn (1984b) c0
    inline constexpr double C_MAX = 0.65;

    struct BedSourcedFlux {
        double flux;  // suspended bed material transport rate, kg/s
        double rouse; // Z = w_s / kappa u_star, overall u_star (dimensionless)
        double z_prime;  // modified rouse number Z' = w_s / (beta kappa u_star) + phi
        double shields_grain;  // grain-related Shields number theta'
        double c_ref;  // reference concentration, volumetric
        double velocity_scale;  // Manning mean velocity / log-law mean velocity over [z_ref, h]
        bool suspending;  // false if theta' <= theta_c, or if the depth does not exceed z_ref
    };

    inline double settlingVelocity(double d50_mm, double rho_s) {
        // Ferguson and Church (2004) settling velocity estimate
        double gamma = (rho_s - RHO_W) / RHO_W;
        double d50_m = d50_mm / 1000.0;
        double velocity = gamma * g * std::pow(d50_m, 2) / (FC_C1 * NU + std::pow(0.75 * FC_C2 * gamma * g * std::pow(d50_m, 3), 0.5));
        return velocity;
    }

    inline double dStar(double d50_mm, double rho_s) {
        // Nondimensional grain size, van Rijn (1984a) eq. 1
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

        // Implementing van Rijn 1984
        double D_star = dStar(d50_mm, rho_s);
        double shields_crit;  // van Rijn (1984a) Shields fit, Fig. 1
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

    // Effective roughness of a movable plane bed, about 3 D90 (van Rijn 1984a, p. 1447), m
    inline double grainRoughness(double d50_mm) {
        return 3.0 * D90_OVER_D50 * d50_mm / 1000.0;
    }

    // Reference level a = ks with a minimum of 0.01 d (van Rijn 1984b, eq. 37), m
    inline double referenceLevel(double d50_mm, double h) {
        return std::max(grainRoughness(d50_mm), 0.01 * h);
    }

    // Grain-related Shields number from u*' = g^0.5 u / C', C' = 18 log10(12 R_b / 3 D90)
    // (van Rijn 1984a, eqs. 16-18), with the bed hydraulic radius R_b taken as the depth
    inline double grainShields(const hydraulics::HydraulicState& hs, double d50_mm, double rho_s) {
        double ratio = 12.0 * hs.depth / grainRoughness(d50_mm);
        if (ratio <= 1.0) {
            return 0.0;  // flow shallower than the roughness: C' undefined
        }
        double chezy_grain = 18.0 * std::log10(ratio);
        double ustar_grain = std::sqrt(g) * hs.velocity / chezy_grain;
        double gamma = (rho_s - RHO_W) / RHO_W;
        return ustar_grain * ustar_grain / (gamma * g * d50_mm / 1000.0);
    }

    // Reference concentration, van Rijn (1984b) eq. 38, volumetric; T from eq. 2
    inline double vanRijnCref(double d50_mm, double rho_s, double a, double shields_grain, double shields_crit) {
        double D_star = dStar(d50_mm, rho_s);
        double T = (shields_grain - shields_crit) / shields_crit;
        double d50_m = d50_mm / 1000.0;
        double c_ref = 0.015 * (d50_m / a) * (std::pow(T, 1.5) / std::pow(D_star, 0.3));
        return std::min(c_ref, C_MAX);  // physical bound; eq. 38 has none
    }

    // Diffusivity ratio (Schmidt number), van Rijn (1984b) eq. 22, stated for 0.1 < w_s/u* < 1 (clamped)
    inline double schmidtNumber(double ws, double ustar) {
        double r = std::clamp(ws / std::max(ustar, USTAR_MIN), 0.1, 1.0);
        return 1.0 + 2.0 * r * r;
    }

    // Overall correction to the Rouse number, van Rijn (1984b) eq. 34,
    // stated for 0.01 <= w_s/u* <= 1 (clamped)
    inline double phiFactor(double ws, double ustar, double c_ref) {
        double r = std::clamp(ws / std::max(ustar, USTAR_MIN), 0.01, 1.0);
        return 2.5 * std::pow(r, 0.8) * std::pow(c_ref / C_MAX, 0.4);
    }

    inline BedSourcedFlux bedSourcedSuspension(const hydraulics::HydraulicState& hs, double B, double d50_mm, double rho_s) {
        // Suspended sediment flux

        BedSourcedFlux flux{};

        // No grain size, no bed material pathway
        if (d50_mm <= 0.0) {
            return flux;
        }

        // Rouse number for every link, including those below threshold
        double ws = settlingVelocity(d50_mm, rho_s);
        double ustar = std::max(hs.ustar, USTAR_MIN);  // overall bed-shear velocity
        flux.rouse = ws / (KAPPA * ustar);

        // Flow too shallow to hold a suspended layer above the reference level
        double z_ref = referenceLevel(d50_mm, hs.depth);
        if (hs.depth <= z_ref) {
            return flux;
        }

        // Transport stage from the grain-related bed shear, not the total ustar
        double shields_crit = criticalShields(d50_mm, rho_s);
        flux.shields_grain = grainShields(hs, d50_mm, rho_s);
        flux.suspending = flux.shields_grain > shields_crit;
        if (!flux.suspending) {
            return flux;
        }

        double c_ref = vanRijnCref(d50_mm, rho_s, z_ref, flux.shields_grain, shields_crit);

        // Modified rouse number Z' = Z + phi (van Rijn 1984b, eqs. 3, 22, 33, 34)
        double Z = ws / (schmidtNumber(ws, ustar) * KAPPA * ustar);
        double Z_prime = Z + phiFactor(ws, ustar, c_ref);

        // Zero-velocity level of the log law, z0 = 0.033 ks (van Rijn 1984b, eq. 40)
        double z0 = 0.033 * grainRoughness(d50_mm);

        int N = 100; // Vertical discretization
        double h = hs.depth;
        double dz = (h - z_ref) / N;
        double log_integral = 0.0;
        double c_log_integral = 0.0;
        for (int i = 0; i <= N; ++i) {
            double z = z_ref + i * dz;

            // Shape of the log-law velocity profile
            double log_term = std::log(z / z0);

            // Concentration profile, van Rijn (1984b) eq. 19a below mid-depth, 19b above
            double c_rel = (z < 0.5 * h)
                         ? std::pow((std::max(h - z, 0.0) / z) * (z_ref / (h - z_ref)), Z_prime)
                         : std::pow(z_ref / (h - z_ref), Z_prime) * std::exp(-4.0 * Z_prime * (z / h - 0.5));
            double c = rho_s * c_ref * c_rel;

            double w = (i == 0 || i == N) ? 0.5 : 1.0;   // trapezoid weights
            log_integral += w * log_term * dz; // m
            c_log_integral += w * c * log_term * dz; // kg / m^2
        }

        // The grain-roughness log law overestimates the Manning mean velocity, so the
        // profile is scaled to carry the routed discharge over [z_ref, h]: the scaled
        // velocity is u(z) = hs.velocity * (h - z_ref) * log_term / log_integral.
        double u_scale = hs.velocity * (h - z_ref) / log_integral;
        flux.velocity_scale = u_scale * KAPPA / ustar;

        flux.flux = u_scale * c_log_integral * B; // multiply by width to get kg / s
        flux.c_ref = c_ref;
        flux.z_prime = Z_prime;

        return flux;
    }
} // namespace sediment
