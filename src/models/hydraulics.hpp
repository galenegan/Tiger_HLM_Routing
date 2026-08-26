#pragma once
#include <algorithm>
#include <cmath>
#include <utility>
#include <limits>
#include <boost/math/tools/roots.hpp>

namespace hydraulics {

    // Physical constants
    const double g = 9.81;
    const double gamma_sed = 1.65;

    // Lower bounds for stable computations
    const double Q_MIN = 1e-8; // m^3/s
    const double SLOPE_MIN = 1e-5; // m/m
    const double WIDTH_MIN = 0.5; // m

    struct HydraulicState {
        double depth;
        double velocity;
        double ustar;
        double shields;
        double froude;
    };

    inline double solveWidth(double width_a, double width_b, double drain_area_km2) {
        double safe_area = std::max(drain_area_km2, 0.0);
        return std::max(width_a * std::pow(safe_area, width_b), WIDTH_MIN);
    }

    inline double solveDepth(double q, double slope, double n, double B, bool wide_channel) {
        using namespace boost::math::tools;

        double wide_channel_approx = std::pow(n * q / (B * std::sqrt(slope)), 0.6);
        if (wide_channel) {
            return wide_channel_approx;
        } else {
            auto residual = [=](double h) {
                return wide_channel_approx * std::pow(1 + 2 * h / B, 0.4) - h;
            };

            eps_tolerance<double> tol(std::numeric_limits<double>::digits - 2);
            std::uintmax_t max_iter = 10;
            std::pair<double, double> result = toms748_solve(residual, wide_channel_approx, wide_channel_approx * 10.0, tol, max_iter);
            double root = (result.first + result.second) / 2.0;
            return root;
        }
    }

    inline HydraulicState computeHydraulics(double q, double slope, double n, double width_a, double width_b, double drain_area_km2, double d50_mm, bool wide_channel) {
        HydraulicState state;

        double B = solveWidth(width_a, width_b, drain_area_km2);
        double h = solveDepth(q, slope, n, B, wide_channel);
        double v = q / (h * B);
        double ustar = std::sqrt(g * h * slope);
        double froude = v / std::sqrt(g * h);

        double d50_m = d50_mm / 1000.0;
        double shields = std::pow(ustar, 2.0) / (gamma_sed * g * d50_m);

        state.depth = h;
        state.velocity = v;
        state.ustar = ustar;
        state.froude = froude;
        state.shields = shields;

        return state;
    }
} // namespace hydraulics
