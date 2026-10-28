#include "sediment.hpp"

#include <algorithm>

#include "models/hydraulics.hpp"
#include "models/sediment_closures.hpp"

void SedimentResults::resize(size_t n_owned, size_t n_steps_sed) {
    n_sed_steps = n_steps_sed;
    const size_t n = n_owned * n_sed_steps;
    hillslope_flux.assign(n, 0.0f);
    total_flux.assign(n, 0.0f);
    rouse.assign(n, 0.0f);
}

void ComputeSedimentAtLevel(const ModelSetup& setup,
  const Partition& part,
  const RunoffData& supply,
  const std::vector<float>& results,
  size_t n_steps,
  const std::vector<size_t>& nodes_at_level,
  SedimentResults& sed) {

    const ModelConfig& cfg = setup.config;
    const size_t n_sed = sed.n_sed_steps;
    const bool use_supply = setup.supply_info.has_value();

    #pragma omp parallel for schedule(runtime)
    for (size_t i = 0; i < nodes_at_level.size(); ++i) {
        const size_t global = nodes_at_level[i];
        if (!part.owns(global)) continue;
        const size_t local = part.local_of[global];

        const NodeInfo& node = setup.node_map.at(global);
        const SedLinkParams& p = setup.sed_params[global];
        const double A_h = node.params[0];
        const double B = hydraulics::solveWidth(p.width_a, p.width_b, p.drain_area_km2);

        // This link's row in the supply file
        const size_t supply_row = use_supply ? supply.idToIndex.at(node.stream_id) : 0;
        for (size_t s = 0; s < n_sed; ++s) {

            // Hillslope washload first
            double hillslope_sourced = 0.0;  // kg/s, calculated from supply.data in kg m-2 hr-1
            if (use_supply) {
                hillslope_sourced = supply.data[supply_row * supply.nTime + s] * A_h / 3600.0;
            }

            // Adding from the parent links
            for (size_t parent : node.parents) {
                hillslope_sourced += sed.hillslope_flux[part.local_of[parent] * n_sed + s];
            }

            // Resuspension from the bed
            double bed_sourced = 0.0;
            double rouse = 0.0;

            // TODO GCE: double check this...
            if (cfg.bed_source == 1) {
                const size_t step = std::min(static_cast<size_t>(s * cfg.supply_resolution / cfg.dt), n_steps - 1);
                const double q = results[local * n_steps + step];
                const auto hs = hydraulics::computeHydraulics(q, p.channel_slope, p.channel_manning_n, p.width_a, p.width_b,
                    p.drain_area_km2, p.d50_mm, cfg.rho_s, cfg.wide_channel);
                const auto bed = sediment::bedSourcedSuspension(hs, B, p.d50_mm, cfg.rho_s);
                bed_sourced = bed.flux;
                rouse = bed.rouse;
            }

            // Adding together
            const size_t out = local * n_sed + s;
            sed.hillslope_flux[out] = static_cast<float>(hillslope_sourced);
            sed.total_flux[out] = static_cast<float>(hillslope_sourced + bed_sourced);
            sed.rouse[out] = static_cast<float>(rouse);

        }
    }

}