#pragma once

#include <cstddef>
#include <vector>

#include "model_setup.hpp"
#include "partition.hpp"

// Struct to store sediment results for a single chunk
// entry [local * n_sed_steps + s] is link part.global_of[local]
// at supply step s
struct SedimentResults {
    size_t n_sed_steps = 0;
    std::vector<float> hillslope_flux;  // kg/s
    std::vector<float> total_flux;  // kg/s, hillslope + bed sourced
    std::vector<float> rouse;  // Rouse number Z

    void resize(size_t n_owned, size_t n_steps_sed);
};

void ComputeSedimentAtLevel(
  const ModelSetup& setup,
  const Partition& part,
  const RunoffData& supply,
  const std::vector<float>& results,
  size_t n_steps,
  const std::vector<size_t>& nodes_at_level,
  SedimentResults& sed
);