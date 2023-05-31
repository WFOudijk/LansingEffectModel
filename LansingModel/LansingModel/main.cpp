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
				
    if (!p.addBinary) p.strengthOfSelection = 0; // survival probability of binary genes will be equal to 1
    if (!p.addAgeSpecific && !p.addQuality) { // survival probability of age-specific genes will be equal to 1
        p.initAgeSpecificGenes = 1;
    }
    if (p.addInvestmentInRepair) {
        // reset number of gametes per female
        p.numOfGametes = p.maximumAge * p.scalingParameterForNumOfOffspring;
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
        
    for (unsigned t = 0; t < p.tEnd; ++t){
        indVec deadIndividualsVec; // to examine all dead individuals
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals); // mortality round of the adults
        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
        // output
        if (t % p.outputTime == 0) { // to prevent every time step of being outputted
            std::cout << t << "\n";
            //createOutputAgeDeath(t, p, ageAtDeath); // generate data for average death age
            createOutputDeclineInGameteQuality(t, p, deadIndividualsVec);
        }
    }
  
    auto t_now = std::chrono::system_clock::now();
    std::chrono::duration<double> diff_t = t_now - t_start;
    std::cout << "First time loop finished. This took: " << diff_t.count() << " seconds = " << diff_t.count() / 60 << " minutes \n";   
    t_start = t_now;
    std::cout << "Choosing " << p.numOfIndividualsToFollow << " number of males and females to research longitudinal. \n";
    
    // track a certain number of individuals
    pop.setTrackedIndividuals(p, rng);

    // to keep track of the followed individuals. If everyone has died, the simulation can stop
    int numOfIndividualsToFollow = p.numOfIndividualsToFollow * 2;
    
    int count = 0;
    while (numOfIndividualsToFollow > 0) {
        count += 1;
        indVec deadIndividualsVec;
        
        pop.reproduce(p, rng); // reproduce to make offspring
        pop.mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals); // mortality round of the adults

        pop.addOffspring(p, rng); // adding offspring to the adults
        pop.mutationRound(p, rng);
								
        numOfIndividualsToFollow = (int) std::count_if(pop.males.begin(), pop.males.end(), [](auto& m) { return m.tracked==1; });
        numOfIndividualsToFollow += std::count_if(pop.females.begin(), pop.females.end(), [](auto& f) { return f.tracked==1; });
				
    }
    
    pop.simulateAgeAtDeath(p, rng);

    std::ofstream ofs;
    ofs.open("outputAgeAlivePop.txt"); // the output file
    for (size_t i =0 ; i < pop.males.size(); ++i) {
        ofs << pop.males[i].age << " "
        << pop.females[i].age << "\n";
    }
    ofs.close();

    std::cout << "Counter = " << count << "\n";
    // only create output of life expectancy for the remaining individuals
    //createOutputLifeExpectancy(p, pop.males, pop.females);
    // create output for the tracked individuals
    createOutputTrackedIndividuals(p, deadTrackedIndividuals);
				
    // to print the duration of the program to the terminal
    t_now = std::chrono::system_clock::now();
    diff_t = t_now - t_start;
    std::cout << "Finished. The program took: " << diff_t.count() << " seconds = " << diff_t.count() / 60 << " minutes \n";
    t_start = t_now;
    return 0;
}
