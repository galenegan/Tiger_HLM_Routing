#pragma once
#include <vector>
#include <string>

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
// Reads `csv_path` (must be UTF‐8, with comma separators, and a header row as above)
// and returns a vector of SedLinkParams, one per data row.
// Throws std::runtime_error on I/O or parse errors.
// ─────────────────────────────────────────────────────────────────────────────
std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path);