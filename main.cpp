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

void save_results(JsonHandler& json,
                  std::map<std::string, std::string>& parameters,
                  std::string filename,
                  std::string timestamp,
                  std::string average,
                  std::string strength,
                  std::string corr_length,
                  std::string results_path) {
    std::map<std::string, std::string> results;

    results["timestamp"] = timestamp;
    results["potential_average"] = average;
    results["potential_average_units"] = "mm^-2";
    results["potential_strength"] = strength;
    results["potential_strength_units"] = "no_units";

    if (corr_length != "") {
        results["correlation_length"] = corr_length;
        results["correlation_length_units"] = "mm";
    }

    // Add parameters to result file
    for (const auto& [key, value] : parameters) {
        results[key] = value;
    }

    std::string result_filename = filename + "_" + timestamp;
    json.dump(results, results_path + result_filename + ".json");
    std::cout << json.content << std::endl;

    // Log the run
    std::ofstream log(results_path + "log.txt", std::ios_base::app);
    log << result_filename << std::endl;
    log.close();
}

int main(void) {
    JsonHandler json;
    std::string timestamp = std::to_string(std::time(0));

    json.load("parameters.json");
    std::string filename = json.json["potential_filename"];
    std::string potential_file = json.json["potential_path"] + filename + CSV;
    u_l_long seed = std::stoi(json.json["random_generator_seed"]);
    std::string initial_condition = json.json["initial_condition"];
    std::string results_path = json.json["results_path"];

    int plane_as_initial = initial_condition == "plane";  // 1 if "plane", 0 if else

    BranchedFlow branches(potential_file, seed);
    branches.initialize(plane_as_initial);
    std::cout << "^ Initialized " << filename << std::endl;

    bool ran_rk4 = false, ran_corr = false;

    if (json.json["run_propagation"] == "true") {
        ran_rk4 = branches.rk4_solve();
        std::cout << "\r^ RK4 done" << std::endl;

        branches.save_film(results_path + filename + "_" + timestamp + "_propagation" + "_" + initial_condition + CSV);
        branches.save_scint(results_path + filename + "_" + timestamp + "_scintillation" + "_" + initial_condition + CSV);
    }

    if (json.json["run_correlation"] == "true") {
        int bins = std::stoi(json.json["correlation_bins"]);
        ran_corr = branches.corr_solve(bins);
        std::cout << "\r^ Corr. done" << std::endl;

        branches.save_corr(results_path + filename + "_correlation" + CSV);
    }

    branches.save_potential(results_path + filename + "_potential"+ CSV);

    save_results(
        json,
        json.json,
        filename,
        timestamp,
        std::to_string(branches.potential_average()),
        std::to_string(branches.potential_strength()),
        (ran_corr) ? std::to_string(branches.correlation_length()) : "",
        results_path);

    return 0;
}