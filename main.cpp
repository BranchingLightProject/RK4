/**
 * Branched Flow RK4 main function
 *
 * All initial conditions are handled by a `parameters.json` file
 */
#include <ctime>
#include <iostream>
#include <map>
#include <string>

#include "deps/branched_flow.h"
#include "deps/json_handler.h"

std::string CSV = ".csv";

int main(void) {
    JsonHandler json;
    std::string timestamp = std::to_string(std::time(0));

    // Load parameters
    json.load("parameters.json");
    std::string filename = json.json["potential_filename"];
    std::string potential_file = json.json["potential_path"] + filename + CSV;
    u_l_long seed = std::stoi(json.json["random_generator_seed"]);
    std::string initial_condition = json.json["initial_condition"];
    std::string results_path = json.json["results_path"];

    BranchedFlow branches(potential_file, seed);
    int plane_as_initial = initial_condition == "beam";  // 0 if "beam", 1 if else
    branches.initialize(plane_as_initial);
    std::cout << "^ Initialized " << filename << std::endl;

    bool ran_rk4 = false, ran_corr = false;

    if (json.json["run_propagation"] == "true") {
        ran_rk4 = branches.rk4_solve();
        std::cout << "\r^ RK4 done" << std::endl;

        branches.save_film(results_path + timestamp + "_propagation" + CSV);
        branches.save_scint(results_path + timestamp + "_scintillation" + CSV);
    }

    if (json.json["run_correlation"] == "true") {
        int bins = std::stoi(json.json["correlation_bins"]);
        ran_corr = branches.corr_solve(bins);
        std::cout << "\r^ Corr. done" << std::endl;

        branches.save_corr(results_path + filename + "_correlation" + CSV);
    }

    branches.save_potential(results_path + filename + CSV);

    std::map<std::string, std::string> results;

    results["timestamp"] = timestamp;
    results["potential_average"] = std::to_string(branches.potential_average());
    results["potential_average_units"] = "mm^-2";
    results["potential_strength"] = std::to_string(branches.potential_strength());
    results["potential_strength_units"] = "no_units";

    if (ran_corr) {
        results["correlation_length"] = branches.correlation_length();
        results["correlation_length_units"] = "mm";
    }

    // Add parameters to result file
    for (const auto& [key, value] : json.json) {
        results[key] = value;
    }

    json.dump(results, results_path + timestamp + ".json");
    std::cout << json.content << std::endl;

    return 0;
}