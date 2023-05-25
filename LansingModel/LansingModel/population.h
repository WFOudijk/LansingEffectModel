//
//  population.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once
#include <iostream>
//#include <oneapi/dpl/algorithm>
//#include <oneapi/dpl/execution> // for parallelization
#include <execution>
#include "individual.h"

using indVec = std::vector<Individual>;

struct Population{
    indVec males;
    indVec females;
    indVec offspring;
    
    void makePopulation(const Parameters& p, Randomizer &rng);
    void reproduce(const Parameters& p, Randomizer& rng);
    void mortalityRound(const Parameters& p,
                        Randomizer& rng,
                        indVec& ageAtDeath,
                        indVec& trackedIndividuals);
    void addOffspring(const Parameters& p, Randomizer& rng);
    void mutationRound(const Parameters& p, Randomizer& rng);
    void setTrackedIndividuals(const Parameters& p, Randomizer& rng);
    void simulateAgeAtDeath(Parameters& p, Randomizer& rng);
    void mortalityRoundOffspring(const Parameters& p, Randomizer& rng, indVec& deadIndividuals);
};

void Population::makePopulation(const Parameters& p,
                                Randomizer &rng){
    /**This function initialises the population and initialises males and females. **/
				
    females.reserve(p.populationSize);
    males.reserve(p.populationSize);
    offspring.reserve(p.populationSize * (p.numOfOffspringPerFemale + 1));
				
    for (size_t i = 0; i < p.populationSize; ++i) {
        males.emplace_back(p, rng, false);
        females.emplace_back(p, rng, true);
    }
}

void Population::reproduce(const Parameters& p,
                           Randomizer& rng){
    /**This function is the reproducing step of the population.  Every female reproduces a numOfOffspringPerFemale
     number of offspring with random males. **/
				
    offspring.clear(); // to make sure the vector is empty
        
    for (auto& f : females){ // loop through every female
        unsigned numOfOffspringPerFemale = p.numOfOffspringPerFemale; // default number of offspring per female

        // checks if investment into repair/ reproduction is included in the model
        if (p.addInvestmentInRepair) numOfOffspringPerFemale = f.calcNumberOfOffspring(p, rng);

        // choose the male to mate with
        auto& mate = males[rng.drawRandomNumber(males.size())];

        // start loop to generate offspring
        for (unsigned i = 0; i < numOfOffspringPerFemale; ++i){ // loop through number of offspring to produce
            offspring.emplace_back(f, mate, rng, p);

            // keep track of offspring of the flagged individuals
            if (mate.tracked) mate.offspring.push_back(offspring.back());
            if (f.tracked) f.offspring.push_back(offspring.back());
        }
    }
}

void Population::mortalityRound(const Parameters& p,
                                Randomizer& rng,                                
                                indVec& deadIndividualsVec,
                                indVec& trackedIndividuals){
    /**This function kills off adults.**/
				
    std::for_each(std::execution::par, begin(females), end(females),
                  [&](auto& f){f.dies(rng,p);});
    std::for_each(std::execution::par, begin(males), end(males),
                  [&](auto& m){ m.dies(rng,p); });
    
    for (size_t m = 0; m < males.size();){
        if (males[m].isDead){ // if this is the case, remove the male from the vector
            if (males[m].tracked) { // the individual is tracked
                trackedIndividuals.push_back(males[m]); // TODO: make move
            }
            deadIndividualsVec.emplace_back(std::move(males[m]));
            males[m] = std::move(males.back());
            males.pop_back();
        } else { // else, continue loop
            ++m;
        }
    }
   
    // same for the females
    for (size_t f = 0; f < females.size();){
        if (females[f].isDead){ // if this is the case, remove female from vector
            if (females[f].tracked){
                trackedIndividuals.push_back(females[f]);
            }
            deadIndividualsVec.emplace_back(std::move(females[f]));
            females[f] = std::move(females.back());
            females.pop_back();
        } else { // else, continue loop
            ++f;
        }
    }
}

void Population::addOffspring(const Parameters& p,
                              Randomizer& rng){
    /**This function adds (random) offspring to the adult vectors untill the vectors are at their maximum size again. **/
				
    while (males.size() < (p.populationSize)){
        int randIndex = rng.drawRandomNumber(offspring.size());
        males.emplace_back(std::move(offspring[randIndex])); // add a random offspring to the males vector
        males.back().makeStemcells(p); // new male has to make stem cells
        males.back().isFemaleSex = 0; // set sex bool to false to indicate this new individual is male
								
        // make sure to remove the offspring to prevent repetition
        offspring[randIndex] = std::move(offspring.back());
        offspring.pop_back();
    }
    
    // same for females
    while (females.size() < (p.populationSize)){
        int randIndex = rng.drawRandomNumber(offspring.size());
        females.emplace_back(std::move(offspring[randIndex])); // add a random offspring to the females vector
        females.back().makeSeveralGametes(p, rng); // new female has to make gametes
        females.back().isFemaleSex = 1; // set sex bool to true to indicate new individual is female
								
        // make sure to remove the offspring to prevent repetition
        offspring[randIndex] = std::move(offspring.back());
        offspring.pop_back();
    }
}

void Population::mutationRound(const Parameters& p,
                               Randomizer &rng){
    /**Function to mutate the gametes from the females and the stem cells from the males. */
				
    std::for_each(std::execution::par,begin(females),end(females),[&](auto& f){ f.mutateGametes(p,rng);});
    std::for_each(std::execution::par,begin(males),end(males),[&](auto& m){ m.mutateGametes(p,rng);});
}

void Population::setTrackedIndividuals(const Parameters &p, Randomizer &rng){
    /**Function to shuffle the remaing population and track a certain number of individuals longitudinally. **/
				
    // shuffle to make it random
    std::shuffle(males.begin(), males.end(), rng.rng);
    std::shuffle(females.begin(), females.end(), rng.rng);
    for (unsigned i = 0; i < p.numOfIndividualsToFollow; ++i){
        // Choosing individuals in order after shuffling the vector to prevent drawing random indviduals without replacement
        males[i].tracked = 1;
        females[i].tracked = 1;
    }
}

void Population::simulateAgeAtDeath(Parameters& p, Randomizer& rng){
    /**Function to simulate age at death from offspring to determine offspring lifespan. **/
    
    // reset the number of offspring per female
    p.numOfOffspringPerFemale = 10;
    
    // female needs more gametes to be able to get more offspring
    for (auto &female : females){
        // every female needs to have enough gametes to produce 10 offspring to track.
        while(female.gametes.size() < p.numOfOffspringPerFemale){
            // copy the gametes into the same vector.
            female.gametes.insert(female.gametes.end(), female.gametes.begin(), female.gametes.end());
        }
    }
    
    // reserve some more memory for the offspring vector
    offspring.reserve(p.populationSize * (p.numOfOffspringPerFemale + 1));

    // have the individuals make offspring
    reproduce(p, rng);
    
    std::ofstream ofs;
    ofs.open("testFile.txt");
    ofs << "Reproduce is finished. \n";
    ofs.close();
    // make new vector to keep track of the dead individuals
    indVec deadIndividuals;
    // reserve memory space for the dead individuals
    deadIndividuals.reserve(offspring.size());
    
    // have offspring go through mortality until they are all dead
    while (!offspring.empty()) {
        // offspring go through mortality round
        mortalityRoundOffspring(p, rng, deadIndividuals);
    }
    
    std::ofstream ofs2;
    ofs2.open("testFile2.txt");
    ofs2 << "mortality is finished. \n";
    
    // make output of the dead individuals
    outputForSimulatedLifespan(deadIndividuals);
    
    std::ofstream ofs3;
    ofs3.open("testFile3.txt");
    ofs3 << "output is finished. \n";
    
}

void Population::mortalityRoundOffspring(const Parameters& p,
                                         Randomizer& rng,
                                         indVec& deadIndividuals){
    /**Function to have the simulated offspring go through a mortality round. **/
    
    std::for_each(std::execution::par, begin(offspring), end(offspring),
                  [&](auto& ind){ind.dies(rng,p);});
       
    // loop through the offspring
    for (size_t indiv = 0; indiv < offspring.size();){
        if (offspring[indiv].isDead){ // if the individual is dead, remove the individual from the vector
            deadIndividuals.emplace_back(std::move(offspring[indiv]));
            offspring[indiv] = std::move(offspring.back());
            offspring.pop_back();
        } else { // else, continue loop
            ++indiv;
        }
    }
}
