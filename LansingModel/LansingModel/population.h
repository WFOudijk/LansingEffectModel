//
//  population.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once
#include <iostream>
#include <optional>
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
};

void Population::makePopulation(const Parameters& p,
                                Randomizer &rng){
    /**This function initialises the population and initialises males and females. **/
    females.reserve(p.populationSize);
    males.reserve(p.populationSize);
    offspring.reserve(p.populationSize * p.numOfOffspringPerFemale);
    for (size_t i = 0; i < p.populationSize; ++i){
	    males.emplace_back(p, rng, false);
        females.emplace_back(p, rng, true);
    }
}

void Population::reproduce(const Parameters& p,
                           Randomizer& rng){
    /**This function is the reproducing step of the population.  Every female reproduces a numOfOffspringPerFemale
     number of offspring with random males. **/
    offspring.clear(); // to make sure the vector is empty
    
    for (auto j = 0u; j < females.size(); ++j){ // loop through every female
        for (unsigned i = 0; i < p.numOfOffspringPerFemale; ++i){ // loop through number of offspring to produce
            int male = rng.drawRandomNumber(males.size());
            Individual newOffspring = Individual(females[j], males[male], rng, p);
            //offspring.emplace_back(females[j], males[rng.drawRandomNumber(males.size())], rng, p); // reproduce
            offspring.push_back(newOffspring);
            if (males[male].identifier) males[male].offspringOfIndividual.push_back(newOffspring);
            if (females[j].identifier) females[j].offspringOfIndividual.push_back(newOffspring);
        }
    }
}

void Population::mortalityRound(const Parameters& p,
                                Randomizer& rng,                                
                                indVec& deadIndividualsVec,
                                indVec& trackedIndividuals){
    /**This function kills off adults.**/
    for (size_t male = 0; male < males.size();){
        bool die = males[male].dies(rng, p); // check if current male will die
        if (die){ // if this is the case, remove the male from the vector
            if (males[male].identifier) {
                males[male].sex = 'M';
                trackedIndividuals.push_back(males[male]);
            }
            deadIndividualsVec.push_back(males[male]);
            males[male] = std::move(males.back());
            males.pop_back();
        } else { // else, continue loop
            ++male;
        }
    }
   
    // same for the females
    for (size_t female = 0; female < females.size();){
        bool die = females[female].dies(rng, p); // check if current female will die
        if (die){ // if this is the case, remove female from vector
            //deadIndividualsVec.push_back(females[female]);
            if (females[female].identifier){
                females[female].sex = 'F';
                trackedIndividuals.push_back(females[female]);
            }
            deadIndividualsVec.push_back(females[female]);
            females[female] = std::move(females.back()); 
            females.pop_back();
        } else { // else, continue loop
            ++female;
        }
    }
}

void Population::addOffspring(const Parameters& p,
                              Randomizer& rng){
    /**This function adds (random) offspring to the adult vectors untill the vectors are at their maximum size again. **/
    while (males.size() < (p.populationSize)){
        int randIndex = rng.drawRandomNumber(offspring.size());
        offspring[randIndex].makeStemcells(p, rng);
        males.push_back(offspring[randIndex]); // add a random offspring to the males vector
        offspring[randIndex] = std::move(offspring.back());
        offspring.pop_back(); // make sure to remove the offspring to prevent repetition
    }
    
    // same for females
    while (females.size() < (p.populationSize)){
        int randIndex = rng.drawRandomNumber(offspring.size());
        offspring[randIndex].makeSeveralGametes(p, rng);
        females.push_back(offspring[randIndex]); // add a random offspring to the females vector
        offspring[randIndex] = std::move(offspring.back());
        offspring.pop_back(); // make sure to remove the offspring to prevent repetition
    }
}

void Population::mutationRound(const Parameters& p,
                               Randomizer &rng){
    for (size_t i = 0; i < females.size(); ++i){
        females[i].mutateGametes(p, rng);
        males[i].mutateStemCells(p, rng); // females and males need to remain same size
     

    }
}

