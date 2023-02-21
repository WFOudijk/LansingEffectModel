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

using arrayOfGenes = std::array<bool, numOfGenes>;

struct Individual {
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    double survivalProb; // based on the binary genes
				double quality; // parental quality
    
    // array with genes - binary
    std::array<arrayOfGenes, 2> genetics;
				
				Parameters p;
					
    // array with age specific genes - survival probabilities
				std::array<double, 40> genesMaternal; // TODO: make this also array of two arrays
				std::array<double, 40> genesPaternal; // TODO: how to get member variable of parameters instead of 40 hardcoded 
    
    // averaging the maternal and paternal survival probabilities
				std::array<double, 40> averageSurvivalProb;
    
    std::vector<Gamete> gametesOfIndividual;
    std::vector<std::array<Gamete, 2> > stemCells;
    
    bool identifier;
    std::vector<Individual> offspringOfIndividual;
    char sex;
    
    Individual(const Parameters& p, // initializing constructor
               Randomizer& rng,
               bool isFemale) : age(0),
                                ageOfMother(0),
                                ageOfFather(0),
																																quality(1.0),
																																identifier(0){
        // initialise two gametes for this individual
								Gamete gameteMaternal(p, rng);
        Gamete gametePaternal(p, rng);
                                    
        // fill the binary gene arrays
        genetics[0] = gameteMaternal.genesOfGamete;
        genetics[1] = gametePaternal.genesOfGamete;
                                    
								// fill the age-specific gene arrays
        genesMaternal = gameteMaternal.ageSpecificGenesOfGamete;
        genesPaternal = gametePaternal.ageSpecificGenesOfGamete;
        
								calcSurvivalProb(p);

        // if the individual is female she should make gametes, otherwise the male should make stem cells.
        (isFemale) ? makeSeveralGametes(p, rng) : makeStemcells(p, rng);
    }

    Individual(Individual& mother,
               Individual& father,
               Randomizer& rng,
               const Parameters& p) : age(0),
                                      identifier(0),
																																				  quality(1.0){
        /**Constructor to reproduce and create offspring . **/
        // first, get a gamete from the mothers gamete list
        Gamete gameteMother = std::move(mother.gametesOfIndividual.back());
        // remove this gamete from the mothers gamete list
        mother.gametesOfIndividual.pop_back(); 
        // have the father generate a gamete
        //Gamete gameteFather = father.makeGamete(rng);
        Gamete gameteFather = father.makeGameteFromStemCell(rng);
        // make a new individual of these gametes
								genetics[0] = gameteMother.genesOfGamete;
								genetics[1] = gameteFather.genesOfGamete;
								genesMaternal = gameteMother.ageSpecificGenesOfGamete;
        genesPaternal = gameteFather.ageSpecificGenesOfGamete;
        // set the ages of the parents
        ageOfMother = mother.age;
        ageOfFather = father.age;
        calcSurvivalProb(p); // to set the survival probability of the new individual
								//survivalProb *= mother.quality; // multiply calculated survival prob with the quality of the mother TODO: needs to be changed.
        }
       
    bool dies(Randomizer& rng, const Parameters& p);
    void calcSurvivalProb(const Parameters& p);
    Gamete makeGamete(Randomizer& rng);
    void makeSeveralGametes(const Parameters &p, Randomizer& rng);
    void mutateGametes(const Parameters& p, Randomizer& rng);
    void makeStemcells(const Parameters& p, Randomizer& rng);
    void mutateStemCells(const Parameters& p, Randomizer& rng);
    Gamete makeGameteFromStemCell(Randomizer& rng);
};

Gamete Individual::makeGamete(Randomizer& rng){
    /**Function to make a single gamete. Based on stochasticity to determine which genes the gamete receives. **/
    Gamete gamete; // initialise gamete
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // based on a random value, determine whether the genes of the gamete come from the mother or from the father
        gamete.genesOfGamete[i] = genetics[rng.bernoulli()][i];
        //gamete.ageSpecificGenesOfGamete[i] = (rng.bernoulli()) ? genesMaternal[i] : genesPaternal[i];
    }
				for (size_t i = 0; i < gamete.ageSpecificGenesOfGamete.size(); ++i){
								gamete.ageSpecificGenesOfGamete[i] = (rng.bernoulli()) ? genesMaternal[i] : genesPaternal[i];
				}
    return gamete;
}

void Individual::makeSeveralGametes(const Parameters &p,
                             Randomizer& rng){
    /** Function to generate multiple gametes for a female.**/
    for (unsigned i = 0; i < p.numOfGametes; ++i){
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
    const double survivalProbForAge = averageSurvivalProb[age] * survivalProb; // the survival prob based on both gene arrays
    double survivalProbIncExtrinsicRisk = survivalProbForAge * (1 - p.extrinsicMortRisk); // taking extrinsic mortality into account
    if (rng.bernoulli(survivalProbIncExtrinsicRisk)){ // bernoulli distribution with the bias of survival probability of the individual
        age += 1; // increment age if individual survives the mortality round
								quality -= p.qualityDecrease; // every time an individual ages, the parental quality should decrease
        if (age == p.maximumAge) dies = true;
    } else { // indidvidual dies
        dies = true; // Individual will die
    }
    return dies;
}

void Individual::mutateGametes(const Parameters &p,
                               Randomizer &rng){
    for (size_t i = 0; i < gametesOfIndividual.size(); ++i){
        gametesOfIndividual[i].mutate(p, rng, false);
    }
    
}
void Individual::makeStemcells(const Parameters& p,
                               Randomizer& rng){
    for (unsigned i = 0; i < p.numOfStemCells; ++i){
        Gamete gamete;
        gamete.genesOfGamete = genetics[0];
								gamete.ageSpecificGenesOfGamete = genesMaternal;
        Gamete gamete2;
        gamete2.genesOfGamete = genetics[1];
								gamete2.ageSpecificGenesOfGamete = genesPaternal;
        std::array<Gamete, 2> genetics = {gamete, gamete2};
        stemCells.push_back(genetics);
    }
}

void Individual::mutateStemCells(const Parameters& p,
                                 Randomizer& rng){
    for (size_t i = 0; i < stemCells.size(); ++i){
        stemCells[i][0].mutate(p, rng, true);
        stemCells[i][1].mutate(p, rng, true);
    }
}

Gamete Individual::makeGameteFromStemCell(Randomizer& rng){
    // first, get a stem cell
    std::array<Gamete, 2> stemCell = std::move(stemCells.back());
    // remove this stem cell from the list
    stemCells.pop_back();
    // initialise gamete
    Gamete gamete;
				// TODO: I use recombination in the generation of the gamete
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // determine based on a bernoulli distribution which gene will be inherited
        gamete.genesOfGamete[i] = stemCell[rng.bernoulli()].genesOfGamete[i];
							//	gamete.ageSpecificGenesOfGamete[i] = stemCell[rng.bernoulli()].ageSpecificGenesOfGamete[i];
    }
				for (size_t i = 0; i < gamete.ageSpecificGenesOfGamete.size(); ++i){
								gamete.ageSpecificGenesOfGamete[i] = stemCell[rng.bernoulli()].ageSpecificGenesOfGamete[i];
				}
    return gamete;
}

void Individual::calcSurvivalProb(const Parameters& p){
				// calculate the survival probability based on the binary represented gene array
				// sum number of ones to calculate the survival probability
				int sumOfDamage1 = std::accumulate(genetics[0].begin(), genetics[0].end(), 0);
				int sumOfDamage2 = std::accumulate(genetics[1].begin(), genetics[1].end(), 0);
				// TODO: determine if both parents quality play a part or at random one or etc
				
				// calculate survival probability based on number of ones
				survivalProb = exp(p.strengthOfSelection * (sumOfDamage1 + sumOfDamage2));
				
				// calculate the array with survival probabilities for the age-specific genes
				for (auto i = 0u; i < genesMaternal.size(); ++i){
								double average = (genesMaternal[i] + genesPaternal[i]) * 0.5;
								averageSurvivalProb[i] = average;
				}
}
