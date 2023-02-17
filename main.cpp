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
#include <iostream>
#include <string>

#include "deps/branched_flow.h"

int INITIAL = 0;
std::string PATH = "../potentials/ising/";
std::string FILENAME = "test";

int main(int argc, char **argv) {
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

    branches.rk4_solve();
    std::cout << "\r^ RK4 done" << std::endl;
    branches.corr_solve(1000);
    std::cout << "\r^ Corr. done" << std::endl;

    branches.save_potential("../results/potentials/" + FILENAME + ".csv");
    branches.save_corr("../results/correlation/" + FILENAME + ".csv");

    if (INITIAL == 1) FILENAME += "_plane";

    branches.save_film("../results/beams/" + FILENAME + ".csv");
    branches.save_scint("../results/scint/" + FILENAME + ".csv");

    std::cout << "\n\n**Results for '" << FILENAME << "'**"
              << "\n- Potential average: " << branches.potential_average() << " [mm^-2]"
              << "\n- Potential strength: " << branches.potential_strength() * 100 << "% []"
              << "\n- Correlation length: " << branches.correlation_length() << " [mm]" << std::endl;

    std::ofstream file("../results/info/" + FILENAME + ".txt");
    file << "\nPotential average: " << branches.potential_average() << " [mm^-2]"
         << "\nPotential strength: " << branches.potential_strength() * 100 << "% []"
         << "\nCorrelation length: " << branches.correlation_length() << " [mm]" << std::endl;
    file.close();

    return 0;
}