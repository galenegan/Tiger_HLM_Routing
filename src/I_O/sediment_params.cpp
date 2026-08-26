#include "sediment_params.hpp"
#include <vector>
#include <string>


std::vector<SedLinkParams> loadSedLinkParams(const std::string& csv_path) {
    std::ifstream in(csv_path);
    if (!in.is_open()) {
        throw std::runtime_error("Failed to open parameter file: " + csv_path);
    }

    std::string line;
    // Read header
    if (!std::getline(in, line)) {
        throw std::runtime_error("Empty parameter file: " + csv_path);
    }
    std::vector<std::string> headers;
    {
        std::istringstream ss(line);
        std::string col;
        while (std::getline(ss, col, ',')) {
            headers.push_back(col);
        }
    }

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
        "drain_area_km2"
        "d50_mm",
    };
    for (auto& name : required) {
        if (idx.find(name) == idx.end()) {
            throw std::runtime_error("Missing column '" + name + "' in " + csv_path);
        }
    }
    std::vector<SedLinkParams> out;
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::vector<std::string> fields;
        std::string cell;
        while (std::getline(ss, cell, ',')) {
            fields.push_back(cell);
        }
        if ((int)fields.size() < (int)headers.size())
            throw std::runtime_error("Bad row with too few fields in " + csv_path);

        SedLinkParams p;
        p.stream_id      = std::stol(fields[idx.at("stream_id")]);
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
        p.drain_area_km2 = std::stod(fields[idx["drain_area_km"]]);
        p.d50_mm = std::stod(fields[idx["d50_mm"]]);

        out.push_back(p);
    }

    return out;
}