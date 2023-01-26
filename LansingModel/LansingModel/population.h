//
//  population.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once
#include <iostream>
#include "individual.h"

using indVec = std::vector<Individual>;

struct Population{
    indVec males;
    indVec females;
    indVec offspring;
    
    void makePopulation(const Parameters& p, Randomizer &rng);
    void reproduce(const Parameters& p, Randomizer& rng);
    void mortalityRound(const Parameters& p, Randomizer& rng, std::vector<Individual>& ageAtDeath);
    void addOffspring(const Parameters& p, Randomizer& rng);
    void mutationRound(const Parameters& p, Randomizer& rng);
};

void Population::makePopulation(const Parameters& p,
                                Randomizer &rng){
    /**This function initialises the population and initialises males and females. **/
    for (int i = 0u; i < p.populationSize; ++i){
	    males.emplace_back(p, rng, false);
        females.emplace_back(p, rng, true);
    }
}

void Population::reproduce(const Parameters& p,
                           Randomizer& rng){
    /**This function is the reproducing step of the population.  Every female reproduces a numOfOffspringPerFemale
     number of offspring with random males. **/
    offspring.clear(); // to make sure the vector is empty
    // to optimize code, reserve the specific space for the offspring vector
    offspring.reserve(females.size() * p.numOfOffspringPerFemale);
    for (auto j = 0u; j < females.size(); ++j){ // loop through every female
        for (int i = 0; i < p.numOfOffspringPerFemale; ++i){ // loop through number of offspring to produce
            offspring.emplace_back(females[j], males[rng.drawRandomNumber(males.size())], rng, p); // reproduce
        }
    }
}

void Population::mortalityRound(const Parameters& p,
                                Randomizer& rng,                                
                                std::vector<Individual>& ageAtDeath){
    /**This function kills off adults.**/
    for (auto male = 0; male < males.size();){
        bool die = males[male].dies(rng, p); // check if current male will die
        if (die){ // if this is the case, remove the male from the vector
            ageAtDeath.push_back(males[male]);
            males[male] = males.back();
            males.pop_back();
        } else { // else, continue loop
            ++male;
        }
    }
   
    // same for the females
    for (auto female = 0; female < females.size();){
        bool die = females[female].dies(rng, p); // check if current female will die
        if (die){ // if this is the case, remove female from vector
            ageAtDeath.push_back(females[female]);
            females[female] = females.back(); // here!
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
        males.push_back(offspring[randIndex]); // add a random offspring to the males vector
        offspring[randIndex] = offspring.back();
        offspring.pop_back(); // make sure to remove the offspring to prevent repetition
    }
    
    // same for females
    while (females.size() < (p.populationSize)){
        int randIndex = rng.drawRandomNumber(offspring.size());
        offspring[randIndex].makeSeveralGametes(p, rng);
        females.push_back(offspring[randIndex]); // add a random offspring to the females vector
        offspring[randIndex] = offspring.back();
        offspring.pop_back(); // make sure to remove the offspring to prevent repetition
    }
}

void Population::mutationRound(const Parameters& p,
                               Randomizer &rng){
    for (auto i = 0; i < females.size(); ++i){
        females[i].mutateGametes(p, rng);
        males[i].mutateStemCells(p, rng); // females and males need to remain same size
     

    }
}

