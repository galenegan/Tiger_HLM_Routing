
#include "model_setup.hpp"  // Your struct and setupModel() declaration

// C++ standard libraries
#include <algorithm>
#include <iostream>         // For std::cout, std::endl
#include <stdexcept>
#include <string>
#include <vector>

// External libraries
#include "I_O/node_info.hpp"
#include "I_O/output_series.hpp"
#include "I_O/inputs.hpp"
#include "I_O/config_loader.hpp"
#include "I_O/sediment_params.hpp"
#include "models/grain_size.hpp"
#include "models/hydraulics.hpp"
#include "models/sediment_closures.hpp"


/**
 * @brief Prints the parsed sediment config, where each link's D50 came from with the spread of
 * the estimated values, and the derived grain properties of the link with node index 0.
 */
static void printSedimentSummary(const ModelSetup& setup) {
    const ModelConfig& c = setup.config;
    std::cout << "  Sediment config:" << std::endl
              << "    wide_channel: " << c.wide_channel
              << ", bed_source: " << c.bed_source
              << ", hillslope_source: " << c.hillslope_source
              << ", rho_s: " << c.rho_s << " kg/m^3" << std::endl
              << "    params_file: " << c.sediment_parameters_file << std::endl
              << "    hydraulics_file: " << c.hydraulics_file << std::endl
              << "    supply_folder: " << c.supply_folder
              << " (varname '" << c.supply_varname
              << "', id_varname '" << c.supply_id_varname
              << "', resolution " << c.supply_resolution << " min)" << std::endl
              << "    sediment_series_filepath: " << c.sediment_series_filepath << std::endl;

    if (setup.sed_params.empty()) return;

    // Where D50 came from, and the spread of the estimated values
    size_t n_provided = 0, n_qbf_column = 0, n_area_law = 0;
    size_t n_sand = 0, n_slope_outside = 0, n_qbf_outside = 0, n_slope_floored = 0;
    std::vector<double> d50_est, p_sand_est;
    for (const SedLinkParams& link : setup.sed_params) {
        if (link.d50_origin == D50Origin::Provided) {
            ++n_provided;
            continue;
        }
        (link.d50_origin == D50Origin::QbfColumn ? n_qbf_column : n_area_law) += 1;
        d50_est.push_back(link.d50_mm);
        p_sand_est.push_back(link.p_sand);
        if (link.p_sand > 0.5) ++n_sand;
        if (link.channel_slope < hydraulics::SLOPE_MIN) ++n_slope_floored;
        const double slope = std::max(link.channel_slope, hydraulics::SLOPE_MIN);
        if (slope < grain_size::SLOPE_CAL_MIN || slope > grain_size::SLOPE_CAL_MAX) ++n_slope_outside;
        if (link.q_bf_m3s < grain_size::QBF_CAL_MIN || link.q_bf_m3s > grain_size::QBF_CAL_MAX) ++n_qbf_outside;
    }
    std::cout << "  D50 source: " << n_provided << " link(s) from d50_mm, " << n_qbf_column
              << " estimated with q_bf_m3s, " << n_area_law << " estimated with Q_bf = "
              << c.q_bf_area_coef << " * A^" << c.q_bf_area_exp << std::endl;
    if (!d50_est.empty()) {
        auto spread = [](std::vector<double> v) {
            std::sort(v.begin(), v.end());
            return std::to_string(v.front()) + " / " + std::to_string(v[v.size() / 2]) + " / " +
                   std::to_string(v.back());
        };
        std::cout << "  Estimated D50 (Tonina et al. 2026), min / median / max: " << spread(d50_est)
                  << " mm; sand-bed probability: " << spread(p_sand_est) << std::endl
                  << "    sand bed (P_s > 0.5): " << n_sand
                  << "; outside calibration range: slope " << n_slope_outside
                  << ", Q_bf " << n_qbf_outside
                  << "; slope floored at " << hydraulics::SLOPE_MIN << ": " << n_slope_floored
                  << std::endl;
    }

    const SedLinkParams& p = setup.sed_params[0];
    std::cout << "  Link at node index 0 (stream_id " << p.stream_id << "):"
              << " d50_mm " << p.d50_mm;
    if (p.d50_mm > 0.0) {
        std::cout << ", settling velocity " << sediment::settlingVelocity(p.d50_mm, c.rho_s) << " m/s"
                  << ", D* " << sediment::dStar(p.d50_mm, c.rho_s)
                  << ", critical Shields " << sediment::criticalShields(p.d50_mm, c.rho_s);
    }
    std::cout << std::endl;
}


/**
 * @brief Sets up the model by loading configuration, node levels, initial conditions,
 * boundary conditions, runoff data, and output options.    
 */
ModelSetup setupModel(const char* config_path) {
    ModelSetup setup;

    // ----------------- SETUP --------------------------------------
    std::cout << "_________________MODEL SET UP_____________________ \n" << std::endl;

    // INPUTS --------------------------------------
    // Read user inputs from YAML file
    std::cout << "Loading user inputs from YAML file...";
    setup.config = ConfigLoader::loadConfig(config_path);
    std::cout << "completed!" << std::endl;
    std::cout << "  snapshot_per_year: " << (setup.config.snapshot_per_year ? "true" : "false") << std::endl;

    // Read node levels from CSV file
    std::cout << "Loading network parameters...";
    read_node_levels(setup.config.parameters_file, setup.node_map, setup.level_groups);
    setup.n_links = setup.node_map.size(); //number of links used for allocating results
    std::cout << "completed!" << std::endl;

    // Load sediment params from CSV file
    std::cout << "Loading sediment parameters...";
    if (setup.config.sediment_flag == 1) {
        const BankfullAreaLaw area_law{setup.config.q_bf_area_coef, setup.config.q_bf_area_exp};
        setup.sed_params = loadSedLinkParams(setup.config.sediment_parameters_file, setup.node_map,
                                             area_law);
        std::cout << "completed!" << std::endl;
        printSedimentSummary(setup);
    } else {
        std::cout << "skipped (sediment.flag 0)" << std::endl;
    }


    // Initial conditions (optional, can be set to a constant value)
    std::cout << "Loading initial conditions...";
    setup.uini = loadInitialConditions(setup.config.initial_conditions_flag,
                                     setup.config.initial_value,
                                     setup.config.initial_conditions_filename,
                                     setup.config.initial_conditions_varname,
                                     setup.config.initial_conditions_id_varname,
                                     setup.n_links);
    std::cout << "completed!" << std::endl;

    //read boundary conditions from file if they exist
    std::cout << "Loading boundary condition...";
    if(setup.config.boundary_conditions_flag == 1){
        setup.boundary_conditions = readBoundaryConditions(setup.config.boundary_conditions_filename,
                                                          setup.config.boundary_conditions_varname,
                                                          setup.config.boundary_conditions_id_varname);
    }
    std::cout << "completed!" << std::endl;

    // Check if reservoir routing is needed
    // This is a placeholder for future implementation
    std::cout << "Checking if reservoir routing is needed (placeholder)...";
    std::cout << "completed!" << std::endl;


    //Get runoff chunk info
    std::cout << "Loading runoff data chunk info...";
    setup.runoff_info = getRunoffChunkInfo(setup.config.runoff_path,
                                           setup.config.runoff_varname,
                                           setup.config.chunk_size);
    std::cout << "completed!" << std::endl;

    // Getting supply chunk info
    if (setup.config.sediment_flag == 1 && setup.config.hillslope_source == 1) {
        std::cout << "Loading supply data chunk info...";
        setup.supply_info = getRunoffChunkInfo(
            setup.config.supply_folder,
            setup.config.supply_varname,
            setup.config.chunk_size
            );

        if (setup.supply_info->nchunks != setup.runoff_info.nchunks) {
            throw std::runtime_error(
                "Supply folder " + setup.config.supply_folder + " gives " +
                std::to_string(setup.supply_info->nchunks) + " chunks but runoff folder " +
                setup.config.runoff_path + " gives " + std::to_string(setup.runoff_info.nchunks));
        }

        for (int i = 0; i < setup.supply_info->nchunks; i++) {
            if (setup.supply_info->ntime[i] != setup.runoff_info.ntime[i]) {
                throw std::runtime_error(
                    "Chunk " + std::to_string(i) + ": supply file " + setup.supply_info->filenames[i] +
                    " gives " + std::to_string(setup.supply_info->ntime[i]) +
                    " time steps but runoff file " + setup.runoff_info.filenames[i] +
                    " gives " + std::to_string(setup.runoff_info.ntime[i]));
            }
        }
        std::cout << "completed!" << std::endl;
    }


    // OUTPUT OPTIONS------------------------
    std::cout << "Setting up output options...";
    if(setup.config.output_flag == 2) setup.save_info = readSaveList(setup.config.link_list_filename); //read the save list from file
    std::cout << "completed!" << std::endl;
    std::cout << "__________________________________________________ \n" << std::endl;

    
    // Print out memory information ----------------------
    // Kills system above 1TB
    size_t input_memory = setup.n_links * setup.config.chunk_size * sizeof(float); // Memory for input runoff data
    size_t results_memory = setup.n_links * setup.config.chunk_size * setup.config.runoff_resolution * sizeof(float) / setup.config.dt; // Memory for results vector
    size_t output_memory = 0.5 * (setup.n_links * setup.config.chunk_size * setup.config.runoff_resolution * sizeof(float)) / setup.config.output_resolution; // Memory for output results
    size_t total_memory = (input_memory + results_memory + output_memory) / (1024.0 * 1024.0 * 1024.0); // Assuming 1.1x for intermediate results and overhead
    if(total_memory > 500.0){
        std::cout << "Warning: Estimated memory usage (saving level 1 and above) is high. Consider reducing chunk size or ensure saving out levels higher than 1 only. Please monitor your system's memory usage." << std::endl;
    }
    return setup;
}