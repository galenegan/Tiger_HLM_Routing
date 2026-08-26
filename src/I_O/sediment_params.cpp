#include "sediment_params.hpp"
#include "node_info.hpp"
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

} // namespace

std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path, std::unordered_map<size_t,NodeInfo>& node_map) {
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
        "d50_mm",
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
        p.channel_manning_n = std::stod(fields[idx["channel_manning_n"]]);
        p.width_a = std::stod(fields[idx["width_a"]]);
        p.width_b = std::stod(fields[idx["width_b"]]);
        p.drain_area_km2 = std::stod(fields[idx["drain_area_km2"]]);
        p.d50_mm = std::stod(fields[idx["d50_mm"]]);

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