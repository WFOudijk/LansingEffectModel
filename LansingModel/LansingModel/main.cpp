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
    //createOutput(pop.females);
    //std::cout << "The males:" << std::endl;
    //createOutput(pop.males);
    
    auto t_start = std::chrono::system_clock::now();
    std::vector<int> ageAtDeath;
    for (int t = 0; t < p.tEnd; ++t){
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, ageAtDeath); // mortality round of the adults
        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
        // output
        //double counter = 0.0;
        //double avgAge = 0.0;
        if (t % p.outputTime == 0) { // to prevent every time step of being outputted
            std::cout << t << std::endl;
//            for (auto i : pop.males) {
//                counter += i.survivalProb;
//                avgAge += i.age;
//            }
            //std::cout << counter / pop.males.size() << std::endl;
            //std::cout << avgAge / pop.males.size() << std::endl;


            //createOutput(pop.males);
            //createOutput(pop.females); 
            createOutputAgeDeath(t, p, ageAtDeath); // generate data for average death age
            //createOuputForGGPlot(pop.males, pop.females, t, p); // generate data for ggplot
        }

    }
    //createOutputLifeExpectancy(pop.males, pop.females, p); // generate data for LE plot

//    createOutput(pop.offspring);

    // to print the duration of the program to the terminal
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_t = t_now - t_start;
    std::cout << "Finished. The program took: " << diff_t.count() << " seconds" << std::endl;
    t_start = t_now;
    return 0;
}
