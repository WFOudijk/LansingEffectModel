//
//  individual.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include "parameters.h"
#include "randomizer.h"
#include "gamete.h"

#include <vector>

struct Individual {
    // the genes represent the survival probabilities;
    std::array<Gamete, 2> genetics;
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    double survivalProb;
    std::vector<Gamete> gametesOfIndividual;
    
    Individual(const Parameters& p,
               Randomizer& rng,
               bool isFemale) : age(0){
        // Initialising constructor. Initialise gene arrays represented by gametes.
        // upon initialisation, the gametes obtain their initial damage.
        Gamete gameteMaternal(p, rng);
        Gamete gametePaternal(p, rng);
        
        genetics[0] = gameteMaternal;
        genetics[1] = gametePaternal;

        calcSurvivalProb(p);

        if (isFemale) makeSeveralGametes(p, rng); // to make a vector of gametes for the females
    }

    Individual(Individual& mother,
               Individual& father,
               Randomizer& rng,
               const Parameters& p) : age(0){
        /**Constructor to reproduce and create offspring . **/
        // first, get a gamete from the mothers gamete list
        Gamete gameteMother = mother.gametesOfIndividual.back();
        // remove this gamete from the mothers gamete list
        mother.gametesOfIndividual.pop_back();
        // have the father generate a gamete
        Gamete gameteFather = father.makeGamete(rng);
        // make a new individual of these gametes
        genetics[0] = gameteMother;
        genetics[1] = gameteFather;
        // set the ages of the parents
        ageOfMother = mother.age;
        ageOfFather = father.age;
        calcSurvivalProb(p); // to set the survival probability of the new individual
        }
       
    bool dies(Randomizer& rng, const Parameters& p);
    void calcSurvivalProb(const Parameters& p);
    Gamete makeGamete(Randomizer& rng);
    void makeSeveralGametes(const Parameters &p, Randomizer& rng);
    void mutateGametes(const Parameters& p, Randomizer& rng);
};

void Individual::calcSurvivalProb(const Parameters& p){
    // sum number of ones to calculate the survival probability
    int sumOfDamage1 = std::accumulate(genetics[0].genesOfGamete.begin(), genetics[0].genesOfGamete.end(), 0);
    int sumOfDamage2 = std::accumulate(genetics[1].genesOfGamete.begin(), genetics[1].genesOfGamete.end(), 0);

    // calculate survival probability based on number of ones
    survivalProb = exp(p.strengthOfSelection * (sumOfDamage1 + sumOfDamage2));
}

Gamete Individual::makeGamete(Randomizer& rng){
    /**Function to make a single gamete. Based on stochasticity to determine which genes the gamete receives. **/
    Gamete gamete; // initialise gamete
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        double pick = rng.uniform(); // first, pick a random value // TODO: make bernoulli 
        // based on this value, determine whether the gene of the gamete is the gene from the mother or from the father
        gamete.genesOfGamete[i] = (pick < 0.5) ? genetics[0].genesOfGamete[i] : genetics[1].genesOfGamete[i];
    }
    return gamete;
}

void Individual::makeSeveralGametes(const Parameters &p,
                             Randomizer& rng){
    /** Function to generate multiple gametes for a female.**/
    for (int i = 0; i < p.numOfGametes; ++i){
        Gamete gamete = makeGamete(rng);
        gametesOfIndividual.push_back(gamete);
    }
}

bool Individual::dies(Randomizer& rng,
                      const Parameters& p){
    /**Function to calculate which individuals will die.**/
    bool dies = false;
    // multiply the survival probability of the individual at its age times
    // the effect of an extrinsic mortality risk on the survival probability
    double survivalProbIncExtrinsicRisk = survivalProb * (1 - p.extrinsicMortRisk);
    if (rng.bernoulli(survivalProbIncExtrinsicRisk)){ // bernoulli distribution with the bias of survival probability of the individual
        age += 1; // increment age if individual survives the mortality round
        if (age == p.maximumAge) dies = true;
    } else { // indidvidual dies
        dies = true; // Individual will die
    }
    return dies;
}

void Individual::mutateGametes(const Parameters &p,
                               Randomizer &rng){
    for (int i = 0; i < gametesOfIndividual.size(); ++i){
        gametesOfIndividual[i].mutate(p, rng);
    }
}
