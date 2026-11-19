#include "sediment_params.hpp"
#include "node_info.hpp"
#include "../models/grain_size.hpp"
#include "../models/hydraulics.hpp"
#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

namespace {

// Strip surrounding whitespace
std::string trimField(const std::string& s) {
    const auto first = s.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return "";
    const auto last = s.find_last_not_of(" \t\r\n");
    return s.substr(first, last - first + 1);
}

// Accept either comma-separated or tab-separated files
char detectDelimiter(const std::string& header_line) {
    return header_line.find('\t') != std::string::npos ? '\t' : ',';
}

std::vector<std::string> splitLine(const std::string& line, char delim) {
    std::vector<std::string> fields;
    std::istringstream ss(line);
    std::string cell;
    while (std::getline(ss, cell, delim)) {
        fields.push_back(trimField(cell));
    }
    return fields;
}

// Value of an optional numeric column; empty when the column is absent or the cell is blank
std::optional<double> optionalField(const std::vector<std::string>& fields,
                                    const std::unordered_map<std::string, int>& idx,
                                    const std::string& name) {
    const auto it = idx.find(name);
    if (it == idx.end() || fields[it->second].empty()) return std::nullopt;
    return std::stod(fields[it->second]);
}

// Fills d50_mm, d50_origin, q_bf_m3s and p_sand by the precedence in sediment_params.hpp
void resolveD50(SedLinkParams& p, std::optional<double> d50_mm, std::optional<double> q_bf_m3s,
                const BankfullAreaLaw& area_law, const std::string& csv_path) {
    const std::string where = "stream_id " + std::to_string(p.stream_id) + " in " + csv_path;
    if (d50_mm) {
        p.d50_mm = *d50_mm;
        p.d50_origin = D50Origin::Provided;
        p.q_bf_m3s = 0.0;
        p.p_sand = -1.0;
        return;
    }
    if (q_bf_m3s) {
        if (!(*q_bf_m3s > 0.0)) {
            throw std::runtime_error(where + ": q_bf_m3s must be positive, got " +
                                     std::to_string(*q_bf_m3s));
        }
        p.q_bf_m3s = *q_bf_m3s;
        p.d50_origin = D50Origin::QbfColumn;
    } else {
        if (!(area_law.coef > 0.0)) {
            throw std::runtime_error(where + " has neither d50_mm nor q_bf_m3s; set "
                                     "sediment.q_bf_area_coef and sediment.q_bf_area_exp to "
                                     "estimate Q_bf from drainage area");
        }
        if (!(p.drain_area_km2 > 0.0)) {
            throw std::runtime_error(where + ": drain_area_km2 must be positive to estimate "
                                     "Q_bf from drainage area, got " +
                                     std::to_string(p.drain_area_km2));
        }
        p.q_bf_m3s = area_law.coef * std::pow(p.drain_area_km2, area_law.exp);
        p.d50_origin = D50Origin::QbfAreaLaw;
    }
    const double slope = std::max(p.channel_slope, hydraulics::SLOPE_MIN);
    const grain_size::BedMaterial bed = grain_size::toninaD50(slope, p.q_bf_m3s);
    p.d50_mm = bed.d50_mm;
    p.p_sand = bed.p_sand;
}

} // namespace

std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path,
                                             const std::unordered_map<size_t, NodeInfo>& node_map,
                                             const BankfullAreaLaw& area_law) {
    std::ifstream in(csv_path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + csv_path);
    }

    std::string line;
    // Read header
    if (!std::getline(in, line)) {
        throw std::runtime_error("Empty parameter file: " + csv_path);
    }
    const char delim = detectDelimiter(line);
    const std::vector<std::string> headers = splitLine(line, delim);

    // Map header name to column index
    std::unordered_map<std::string,int> idx;
    for (int i = 0; i < static_cast<int>(headers.size()); ++i) {
        idx[headers[i]] = i;
    }

    // Required columns: match CSV's names
    const std::vector<std::string> required = {
        "stream_id",
        "channel_slope",
        "channel_manning_n",
        "width_a",
        "width_b",
        "drain_area_km2",
    };
    for (auto& name : required) {
        if (idx.find(name) == idx.end()) {
            // Report what was actually parsed
            std::string found;
            for (const auto& h : headers) found += " '" + h + "'";
            throw std::runtime_error(
                "Missing required column '" + name + "' in " + csv_path +
                ". Columns found:" + found +
                ". Expected a header row with columns separated by commas or tabs.");
        }
    }

    // Building a stream_id -> index map
    std::unordered_map<int, size_t> id_to_idx;
    for (const auto& [key, value] : node_map) {
        id_to_idx[value.stream_id] = key;
    }

    // Initializing output to fill in.
    std::vector<SedLinkParams> out(node_map.size());
    std::vector<bool> filled(out.size(), false);
    for (const auto& [key, value] : node_map) {
        if (key >= out.size()) {
            throw std::runtime_error(
                "node index " + std::to_string(key) + " exceeds the link count (" +
                std::to_string(out.size()) + "); node indices must be contiguous from 0");
        }
    }

    while (std::getline(in, line)) {
        if (trimField(line).empty()) continue;
        const std::vector<std::string> fields = splitLine(line, delim);
        if ((int)fields.size() < (int)headers.size())
            throw std::runtime_error("Bad row with too few fields in " + csv_path);

        SedLinkParams p;
        p.stream_id = std::stol(fields[idx.at("stream_id")]);
        // Optional next_stream: default to -1 when the column is absent or blank
        auto it_next = idx.find("next_stream_id");
        if (it_next != idx.end()) {
            const std::string& v = fields[it_next->second];
            if (!v.empty())
                p.next_stream_id = std::stol(v);
            else
                p.next_stream_id = -1;
        } else {
            p.next_stream_id = -1;
        }

        // direct mappings
        p.channel_slope = std::stod(fields[idx.at("channel_slope")]);
        p.channel_manning_n = std::stod(fields[idx.at("channel_manning_n")]);
        p.width_a = std::stod(fields[idx.at("width_a")]);
        p.width_b = std::stod(fields[idx.at("width_b")]);
        p.drain_area_km2 = std::stod(fields[idx.at("drain_area_km2")]);
        resolveD50(p, optionalField(fields, idx, "d50_mm"), optionalField(fields, idx, "q_bf_m3s"),
                   area_law, csv_path);

        auto it = id_to_idx.find(p.stream_id);
        if (it == id_to_idx.end()) {
            throw std::runtime_error("stream_id " + std::to_string(p.stream_id) +
                               " in " + csv_path + " is not in the routing network");
        }
        const size_t out_index = it->second;
        if (out_index >= out.size()) {
            throw std::runtime_error("node index " + std::to_string(out_index) +
                                     " out of range for " + std::to_string(out.size()) + " links");
        }

        if (filled[out_index]) {
            throw std::runtime_error("duplicate row for stream_id " +
                                     std::to_string(p.stream_id) + " in " + csv_path);
        }
        out[out_index] = p;
        filled[out_index] = true;
    }

    size_t n_missing = 0;
    for (const auto& [index, node] : node_map) {
        if (!filled[index]) {
            ++n_missing;
        }
    }
    if (n_missing > 0) {
        throw std::runtime_error(std::to_string(n_missing) + " network link(s) missing from " + csv_path);
    }

    return out;
}