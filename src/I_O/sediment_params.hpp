#pragma once
#include <vector>
#include <string>
#include <unordered_map>
#include "node_info.hpp"

// Where a link's D50 came from
enum class D50Origin {
    Provided,    // d50_mm column
    QbfColumn,   // Tonina et al. (2026) with Q_bf from the q_bf_m3s column
    QbfAreaLaw,  // Tonina et al. (2026) with Q_bf = coef * A_d^exp
};

struct SedLinkParams {
    int stream_id;
    int next_stream_id;
    double channel_slope;
    double channel_manning_n;
    double width_a;
    double width_b;
    double drain_area_km2;
    double d50_mm;
    D50Origin d50_origin;
    double q_bf_m3s;  // bankfull discharge used to estimate D50, m^3/s; 0 when D50 was provided
    double p_sand;    // sand-bed probability from the estimate; -1 when D50 was provided
};

// Bankfull discharge as a power law of drainage area, Q_bf = coef * A_d^exp (m^3/s, A_d in
// km^2). coef <= 0 means not set.
struct BankfullAreaLaw {
    double coef = 0.0;
    double exp = 0.0;
};

// ─────────────────────────────────────────────────────────────────────────────
// loadSedLinkParams
// ─────────────────────────────────────────────────────────────────────────────
// Reads per-link channel and sediment properties and returns them indexed by
// node.index (NOT by row order), so out[node.index] matches node_map.
//
// File format:
//   - One header row naming the columns, then one row per link.
//   - Fields separated by commas or tabs; the separator is inferred from the
//     header. Either LF or CRLF line endings. Surrounding whitespace is ignored.
//   - Column order does not matter; columns are matched by name.
//   - Required columns:
//       stream_id           routing link ID, must match params.csv
//       channel_slope       bed slope S, m/m
//       channel_manning_n   Manning roughness n, s/m^(1/3)
//       width_a             coefficient of B = a * A_d^b, m per km^(2b)
//       width_b             exponent of B = a * A_d^b, dimensionless
//       drain_area_km2      upstream drainage area A_d, km^2
//   - Optional columns (absent column or blank cell means not provided):
//       next_stream_id      defaults to -1
//       d50_mm              median grain size D50, mm; 0 disables the bed pathway
//       q_bf_m3s            bankfull discharge Q_bf, m^3/s
//
// D50 for each link, in order of precedence:
//   1. d50_mm when provided;
//   2. otherwise grain_size::toninaD50(S, Q_bf) with Q_bf from q_bf_m3s;
//   3. otherwise the same with Q_bf = area_law.coef * A_d^area_law.exp.
// S in the estimate is floored at hydraulics::SLOPE_MIN, as in the hydraulics.
//
// Example (comma-separated):
//   stream_id,channel_slope,channel_manning_n,width_a,width_b,drain_area_km2,d50_mm,q_bf_m3s
//   420558772,0.01,0.035,2.7,0.5,0.115,,0.156
//
// Every link in node_map must appear exactly once. Throws std::runtime_error on
// I/O or parse errors, a missing required column, a stream_id absent from the
// routing network, a duplicate row, any network link missing from the file, a
// provided q_bf_m3s that is not positive, or a link that needs the area law when
// area_law.coef is not set or drain_area_km2 is not positive.
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path,
                                             const std::unordered_map<size_t, NodeInfo>& node_map,
                                             const BankfullAreaLaw& area_law);
