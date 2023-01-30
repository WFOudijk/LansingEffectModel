//
//  main.cpp
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//

#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <string>
#include "individual.h"
#include "parameters.h"
#include "outputGenerator.h"
#include "population.h"

int main(int argc, const char * argv[]) {
    // obtain seed from system clock
    std::chrono::high_resolution_clock::time_point tp =
        std::chrono::high_resolution_clock::now();
    unsigned seed = static_cast<unsigned>(tp.time_since_epoch().count());
    std::cout << "The seed is: " << seed << std::endl;
    // create and seed pseudo-random number generator
    Randomizer rng;
    rng.setSeed(seed);
    
    Parameters p; // make parameters object
    
    // read parameter file
    std::string parameterFile;
    if (argc > 1){
        parameterFile = argv[1];
        p.readParameters(parameterFile); // sets parameters from file
    }

    Population pop;
    pop.makePopulation(p, rng); // initialise population
    
    auto t_start = std::chrono::system_clock::now();
    for (int t = 0; t < p.tEnd; ++t){
        std::vector<Individual> ageAtDeath;
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, ageAtDeath); // mortality round of the adults
        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
        // output
        if (t % p.outputTime == 0) { // to prevent every time step of being outputted
            std::cout << t << std::endl;
            //createOutputAgeDeath(t, p, ageAtDeath); // generate data for average death age
            createOutputDeclineInGameteQuality(t, ageAtDeath);
        }

    }
    // to print the duration of the program to the terminal
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_t = t_now - t_start;
    std::cout << "Finished. The program took: " << diff_t.count() << " seconds" << std::endl;
    t_start = t_now;
    return 0;
}
