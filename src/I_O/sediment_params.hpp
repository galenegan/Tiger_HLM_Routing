#pragma once
#include <vector>
#include <string>
#include "node_info.hpp"

struct SedLinkParams {
    int stream_id;
    int next_stream_id;
    double channel_slope;
    double channel_manning_n;
    double width_a;
    double width_b;
    double drain_area_km2;
    double d50_mm;
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
//       d50_mm              median grain size D50, mm
//   - Optional column: next_stream_id (defaults to -1 when absent or blank).
//
// Example (comma-separated):
//   stream_id,channel_slope,channel_manning_n,width_a,width_b,drain_area_km2,d50_mm
//   420558772,0.01,0.035,2.7,0.5,0.115,0.5
//
// Every link in node_map must appear exactly once. Throws std::runtime_error on
// I/O or parse errors, a missing required column, a stream_id absent from the
// routing network, a duplicate row, or any network link missing from the file.
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path, std::unordered_map<size_t,NodeInfo>& node_map);