//
//  population.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once
#include <iostream>
// disable the following two lines when code needs to run on the cluster
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution> // for parallelization
// enable the following two lines when code needs to run on the cluster
//#include <algorithm>
//#include <execution>
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
    void simulateAgeAtDeath(Parameters& p, Randomizer& rng);
    void mortalityRoundOffspring(const Parameters& p, Randomizer& rng, indVec& deadIndividuals);
    void simulateOffspringLifespan(const Parameters& p, Randomizer& rng);
};

void Population::makePopulation(const Parameters& p,
                                Randomizer &rng){
    /**Initialise the population by initialising males and females. **/
				
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
     number of offspring with a random male. **/
				
    offspring.clear(); // to make sure the vector is empty
        
    for (auto& f : females){ // loop through every female
        
        // choose the male to mate with
        auto& mate = males[rng.drawRandomNumber(males.size())];

        // start loop to generate offspring
        for (unsigned i = 0; i < p.numOfOffspringPerFemale; ++i){ // loop through number of offspring to produce
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
    /**This function kills removes the individuals that have died from the population.**/
				
    std::for_each(std::execution::par, begin(females), end(females),
                  [&](auto& f){f.dies(rng,p);});
    std::for_each(std::execution::par, begin(males), end(males),
                  [&](auto& m){ m.dies(rng,p); });
    
    // loop through every male
    for (size_t m = 0; m < males.size();){
        if (males[m].isDead){ // if male is dead, remove from the vector
            if (males[m].tracked) { // if the individual is tracked ..
                // .. the male information is recorded in this vector ..
                trackedIndividuals.push_back(males[m]);
            }
            // .. otherwise in this vector
            deadIndividualsVec.emplace_back(std::move(males[m]));
            
            // remove male from vector
            males[m] = std::move(males.back());
            males.pop_back();
        } else { // else, continue loop
            ++m;
        }
    }
   
    // same for the females
    for (size_t f = 0; f < females.size();){
        if (females[f].isDead){ // if female is dead, remove from the vector
            if (females[f].tracked){ // if the individual is tracked ..
                // .. the female information is recorded in this vector ..
                trackedIndividuals.push_back(females[f]);
            }
            // .. otherwise in this vector
            deadIndividualsVec.emplace_back(std::move(females[f]));
            
            // remove female from vector
            females[f] = std::move(females.back());
            females.pop_back();
        } else { // else, continue loop
            ++f;
        }
    }
}

void Population::addOffspring(const Parameters& p,
                              Randomizer& rng){
    /**This function adds (random) offspring to the adult vectors for every adult that has died. **/
    
    auto maxCapacityF = p.populationSize; // female max capacity
    auto maxCapacityM = p.populationSize; // male max capacity
    
    // checks if there are not enough offspring to fill both the male and the female vectors to their
    // initial size, then fill the distribute the offspring over the females and males
    if (2 * p.populationSize - (males.size() + females.size()) > offspring.size()){
        auto halfOffspring = offspring.size() * 0.5;
        maxCapacityF = females.size() + halfOffspring; // get new maximum capacity
        maxCapacityM = males.size() + halfOffspring; // same for the males
    }
    
    while (males.size() < maxCapacityM){
        int randIndex = rng.drawRandomNumber(offspring.size());
        males.emplace_back(std::move(offspring[randIndex])); // add a random offspring to the males vector
        males.back().makeStemcells(p); // new male has to make stem cells
        males.back().isFemaleSex = 0; // set sex bool to false to indicate this new individual is male
								
        // make sure to remove the offspring to prevent repetition
        offspring[randIndex] = std::move(offspring.back());
        offspring.pop_back();
    }
    
    // same for females
    while (females.size() < maxCapacityF){
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
    std::for_each(std::execution::par,begin(males),end(males),[&](auto& m){ m.mutateStemCells(p,rng);});
}


void Population::simulateAgeAtDeath(Parameters& p, Randomizer& rng){
    /**Function to simulate age at death from offspring to determine offspring lifespan for the cross-sectional analysis. **/
    
    // save the num of offspring per female in variable
    unsigned int saveNumOfOffspringPerFemale = p.numOfOffspringPerFemale;
    
    // reset the number of offspring per female
    p.numOfOffspringPerFemale = p.numOfOffspringForOffspringLifespanSim;
        
    // reserve some more memory for the offspring vector
    offspring.reserve(p.populationSize * (p.numOfOffspringPerFemale + 1));

    // have the individuals make offspring
    reproduce(p, rng);
    
    // make new vector to keep track of the dead individuals
    indVec deadIndividuals;
    // reserve memory space for the dead individuals
    deadIndividuals.reserve(offspring.size() + 1);
    
    // let offspring live out their lives and record the age at death
    while (!offspring.empty()) {
        // offspring go through mortality round
        mortalityRoundOffspring(p, rng, deadIndividuals);
    }

    // make output of the dead individuals and their offspring
    outputForSimulatedLifespan(deadIndividuals);
    
    // reset the number of offspring per female
    p.numOfOffspringPerFemale = saveNumOfOffspringPerFemale;
    
}

void Population::mortalityRoundOffspring(const Parameters& p,
                                         Randomizer& rng,
                                         indVec& deadIndividuals){
    /**Function to let the offspring live out their lives for the cross-sectional analysis. **/
    
    std::for_each(std::execution::par, begin(offspring), end(offspring),
                  [&](auto& ind){ind.dies(rng,p);});
       
    // loop through the offspring
    for (size_t indiv = 0; indiv < offspring.size();){
        if (offspring[indiv].isDead){ // if the offspring is dead, remove from the vector
            deadIndividuals.emplace_back(std::move(offspring[indiv]));
            offspring[indiv] = std::move(offspring.back());
            offspring.pop_back();
        } else { // else, continue loop
            ++indiv;
        }
    }
}

void Population::simulateOffspringLifespan(const Parameters& p,
                                           Randomizer& rng){
    /**Function to look at the individuals longitudinally to determine offspring lifespan. **/
    
    // population reproduces
    reproduce(p, rng);
    
    auto tmp = offspring.size() * 0.5; // get the half value
    
    // make the offspring the new population
    males = {offspring.begin(), offspring.end() - tmp};  // first half
    females = {offspring.begin() + tmp, offspring.end()}; // second half

    // have males make stem cells
    std::for_each(std::execution::par, begin(males), end(males),
                  [&](auto& ind){
        ind.makeStemcells(p);
        ind.isFemaleSex = 0;
    });

    // have females make gametes
    std::for_each(std::execution::par, begin(females), end(females),
                  [&](auto& ind){
        ind.makeSeveralGametes(p, rng);
        ind.isFemaleSex = 1;
        // set tracked to true to keep track of the offspring from the females
        ind.tracked = 1;
    });
    
    // needed for the mortality function - not used in this case
    indVec deadIndividualsVec;
    // to keep track of the females and their offspring over their lives
    indVec deadTrackedIndividuals;
    
    // loop until maximum age.
    for (size_t i = 0; i < p.maximumAge; ++i){
        
        // only reproduce if they are both not empty
        if (!males.empty() && !females.empty()) reproduce(p, rng); // reproduce to make offspring
        
        // mortality round of the adults
        mortalityRound(p, rng, deadIndividualsVec, deadTrackedIndividuals);
        
        // mutation round of the gametes/ stem cells
        mutationRound(p, rng);
    }
    
    // loop through the offspring of the tracked individuals and let them live out their life
    for (Individual& ind : deadTrackedIndividuals){
        // use the counter to track which of the offspring are already dead
        size_t counter = 0;
        // while the counter is not equal to the offspring vector size, there are still offspring alive.
        while (counter != ind.offspring.size()){
            counter = 0;
            for (Individual& offs : ind.offspring){
                // checks if the offspring is already dead
                if (!offs.isDead){
                    // if not, go through mortality round, it either dies or ages one year.
                    offs.dies(rng, p);  
                } else {
                    counter += 1;
                }
            }
        }
    }
    
    // create output for offspring lifespan for longitudinal analysis
    outputOffspringLifespanLongitudinal(deadTrackedIndividuals);
    
    // create output for with age-specific gene values
    createOutputAgeSpecificGenes(p, deadTrackedIndividuals);
}
