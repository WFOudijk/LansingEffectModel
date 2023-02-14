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

    // set the mutationEffect distribution with mean and sd of mutation
    rng.setMutationEffect(p.meanMutationBias, p.sdMutationalEffectSize);
    
    Population pop;
    pop.makePopulation(p, rng); // initialise population
    
    indVec deadTrackedIndividuals; // TODO: should be a better solution

    auto t_start = std::chrono::system_clock::now();
    for (int t = 0; t < p.tEnd; ++t){
        std::vector<Individual> deadIndividualsVec; // to examine all dead individuals 
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals); // mortality round of the adults
        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
        // output
        if (t % p.outputTime == 0) { // to prevent every time step of being outputted
            std::cout << t << std::endl;
            //createOutputAgeDeath(t, p, ageAtDeath); // generate data for average death age
            createOutputDeclineInGameteQuality(t, p, deadIndividualsVec);
        }

    }
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_t = t_now - t_start;
    std::cout << "First time loop finished. This took: " << diff_t.count() << " seconds = " << diff_t.count() / 60 << " minutes " << std::endl;
    t_start = t_now;
    std::cout << "Choosing " << p.numOfIndividualsToFollow << " number of males and females to research longitudinal." << std::endl;
    std::default_random_engine rand(seed);
        
    // shuffle to make it random
    std::shuffle(pop.males.begin(), pop.males.end(), rand);
    std::shuffle(pop.females.begin(), pop.females.end(), rand);
    for (unsigned i = 0; i < p.numOfIndividualsToFollow; ++i){
        // Choosing individuals in order after shuffling the vector to prevent drawing random indviduals without replacement
        pop.males[i].identifier = 1;
        pop.females[i].identifier = 1;
    }
    // to keep track of the followed individuals. If everyone has died, the simulation can stop
    int numOfIndividualsToFollow = p.numOfIndividualsToFollow * 2;
    
    while (numOfIndividualsToFollow > 0) {
        indVec deadIndividualsVec;
        
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals); // mortality round of the adults

        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
        
        numOfIndividualsToFollow = 0;
        for (size_t ind = 0; ind < pop.males.size(); ++ind){
            if (pop.males[ind].identifier) numOfIndividualsToFollow += 1;
            if (pop.females[ind].identifier) numOfIndividualsToFollow += 1;
        }
    }
        
    // only create output of life expectancy for the remaining individuals
    createOutputLifeExpectancy(p, pop.males, pop.females);
    // create output for the tracked individuals
    createOutputTrackedIndividuals(p, deadTrackedIndividuals);
    // to print the duration of the program to the terminal
    t_now = std::chrono::system_clock::now();
    diff_t = t_now - t_start;
    std::cout << "Finished. The program took: " << diff_t.count() << " seconds = " << diff_t.count() / 60 << " minutes " << std::endl;
    t_start = t_now;
    return 0;
}
