//
//  main.cpp
//  Parental age and offspring lifespan: the Lansing effect and its underlying mechanisms.
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//  s4995805
//  Master Biology: Modelling in the life sciences.

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
    std::cout << "The seed is: " << seed << "\n" ;
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
    
    p.setAdditionalParams();
    
    // disable the effects of these mechanisms if they are off in the model run
    if (!p.addBinary) p.strengthOfSelection = 0; // survival probability of binary genes will be equal to 1
    if (!p.addAgeSpecific && !p.addQuality) { // survival probability of age-specific genes will be equal to 1
        p.initAgeSpecificGenes = 1;
    }    

    // set the mutationEffect distribution with mean and sd of mutation
    rng.setMutationEffect(p.meanMutationBias, p.sdMutationalEffectSize);
    // set the number of events distribution for the age specific genes
    rng.setDistMutEvents(p.mutationProbAgeSpecificGenes * p.maximumAge);
    // set the mutation effect for age-specific genes involved in repair/reproduction
    rng.setMutationEffectInvestment(p.meanMutationBiasInvestmentInRepair, p.sdMutationalEffectInvestmentInRepair);
    
    Population pop;
    pop.makePopulation(p, rng); // initialise population
    
    indVec deadTrackedIndividuals;

    auto t_start = std::chrono::system_clock::now();
    
    // start simulation
    for (unsigned t = 0; t < p.tEnd; ++t){
        indVec deadIndividualsVec; // to examine all dead individuals
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals); // mortality round of the adults
        pop.addOffspring(p, rng); // adding offspring to replace dead adults
        pop.mutationRound(p, rng); // mutate the gametes and stem cells
        // output
        if (t % p.outputTime == 0) { // to prevent every time step of being outputted
            std::cout << t << "\n";
            //createOutputDeclineInGameteQuality(t, p, deadIndividualsVec);
        }
    }
  
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_t = t_now - t_start;
    std::cout << "Time loop finished. This took: " << diff_t.count() << " seconds = " << diff_t.count() / 60 << " minutes \n";
    
    // to record the ages of the population after the simulation has ended
    std::ofstream ofs;
    ofs.open("outputAgeAlivePop.txt"); // the output file
    for (size_t i =0 ; i < pop.males.size(); ++i) {
        ofs << pop.males[i].age << " "
        << pop.females[i].age << "\n";
    }
    ofs.close();
    
    // simulate cross-sectional offspring lifespan
    pop.simulateAgeAtDeath(p, rng);
    
    // simulate longitudinal offspring lifespan over parental age
    pop.simulateOffspringLifespan(p, rng);
    				
    return 0;
}
