#pragma once
#include <cmath>

// Median bed-material grain size from channel slope and bankfull discharge,
// Tonina, Villarini and Marzadri (2026), "Distinct Sand- and Gravel-Bed Regimes Structure
// Riverbed Grain-Size Scaling", eqs. 31-34 and Table 1.
// Units: S dimensionless, Q_bf in m^3/s, D50 in mm.

namespace grain_size {

    // Calibration ranges (Sec. 2.1); outside them the relation is an extrapolation (Sec. 2.2.3)
    inline constexpr double SLOPE_CAL_MIN = 8.75e-6;
    inline constexpr double SLOPE_CAL_MAX = 0.0959;
    inline constexpr double QBF_CAL_MIN = 0.021;   // m^3/s
    inline constexpr double QBF_CAL_MAX = 2.16e5;  // m^3/s

    struct BedMaterial {
        double d50_mm;  // median bed-material grain size, mm
        double p_sand;  // probability of a sand bed, P(D50 < 2 mm)
    };

    // Sand-bed probability, eq. 31: logit(P_s) = ln(P_s / (1 - P_s)) (eq. 2)
    inline double sandBedProbability(double slope, double q_bf) {
        double logit = -12.1362 - 4.2675 * std::log10(slope) - 0.6931 * std::log10(q_bf);
        return 1.0 / (1.0 + std::exp(-logit));
    }

    // D50 = D50_s^P_s * D50_g^(1 - P_s), eq. 34, with the slope-only sand branch (eq. 32) and the
    // slope-discharge gravel/coarser branch (eq. 33). Requires slope > 0 and q_bf > 0.
    inline BedMaterial toninaD50(double slope, double q_bf) {
        double p_sand = sandBedProbability(slope, q_bf);
        double d50_sand = 0.9291 * std::pow(slope, 0.1295);
        double d50_gravel = 258.47 * std::pow(slope, 0.5212) * std::pow(q_bf, 0.3058);
        return {std::pow(d50_sand, p_sand) * std::pow(d50_gravel, 1.0 - p_sand), p_sand};
    }

} // namespace grain_size
