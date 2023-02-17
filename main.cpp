/**
 * Branched Flow RK4 main function
 *
 * If you pass command line arguments to this program it will override the path and filename (and perhaps
 * the initial wave profile too) defined in here.
 *
 * - The second argument (the first if you don't include the program name) will be the path
 * - The third argument will be the filename
 * - If you pass a fourth argument (it doesn't matter which one) the initial wave profile will be a plane wave
 */
#include <ctime>
#include <iostream>
#include <map>
#include <string>

#include "deps/branched_flow.h"
#include "deps/json_handler.h"

int INITIAL = 0;
std::string PATH = "../potentials/ising/";
std::string FILENAME = "test";
std::string RESULTS = "./";

int main(int argc, char **argv) {
    JsonHandler json;
    std::string timestamp = std::to_string(std::time(0));

    // TODO: initial conditions & run controls using json

    if (argc > 1) {
        PATH = argv[1];
        FILENAME = argv[2];

        if (argc > 3)
            INITIAL = 1;
        else
            INITIAL = 0;
    }

#if POTENTIAL_SOURCE == 0
    BranchedFlow branches(PATH + FILENAME + ".csv", 1);
    branches.initialize(INITIAL);
    std::cout << "^ Initialized from file " << FILENAME << std::endl;
#else
    BranchedFlow branches(1);
    branches.initialize(INITIAL);
    std::cout << "^ Initialized from function" << std::endl;
#endif

    bool ran_rk4 = false, ran_corr = false;

    ran_rk4 = branches.rk4_solve();
    std::cout << "\r^ RK4 done" << std::endl;
    ran_corr = branches.corr_solve(1000);
    std::cout << "\r^ Corr. done" << std::endl;

    branches.save_potential("../results/potentials/" + FILENAME + ".csv");
    branches.save_corr("../results/correlation/" + FILENAME + ".csv");
    branches.save_film("../results/beams/" + FILENAME + ".csv");
    branches.save_scint("../results/scint/" + FILENAME + ".csv");

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

    json.dump(results, RESULTS + timestamp + ".json");
    std::cout << json.content << std::endl;

    return 0;
}