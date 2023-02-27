//
//  individual.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include "gamete.h"

#include <vector>

using arrayOfGenes = std::array<bool, numOfGenes>;

struct Individual {
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    float survivalProb; // based on the binary genes
				float quality; // parental quality
				bool tracked; // the flagged individuals will be true
				
    // array with genes - binary
    std::array<arrayOfGenes, 2> geneticsBinary;
				
				static Parameters p;
					
    // array with age specific genes - survival probabilities
				std::array<float, 40> ageSpecificGenesMaternal; // TODO: make this also array of two arrays
				// use vector instead of array and reserve space

				std::array<float, 40> ageSpecificGenesPaternal;
    
    // averaging the maternal and paternal survival probabilities
				std::array<float, 40> averageSurvivalProbAgeGenes;
    
    std::vector<Gamete> gametes;
    std::vector<std::array<Gamete, 2> > stemCells;
    
    std::vector<Individual> offspringOfIndividual;
    bool isFemaleSex;
    
    Individual(const Parameters& p, // initializing constructor
               Randomizer& rng,
															bool isFemale);

    Individual(Individual& mother,
               Individual& father,
               Randomizer& rng,
												   const Parameters& p);
       
    bool dies(Randomizer& rng, const Parameters& p);
    void calcSurvivalProb(const Parameters& p);
    Gamete makeGamete(Randomizer& rng);
    void makeSeveralGametes(const Parameters &p, Randomizer& rng);
    void mutateGametes(const Parameters& p, Randomizer& rng);
    void makeStemcells(const Parameters& p);
    void mutateStemCells(const Parameters& p, Randomizer& rng);
    Gamete makeGameteFromStemCell(Randomizer& rng);
};

Individual::Individual(const Parameters& p, // initializing constructor
																							Randomizer& rng,
																							bool isFemale) : age(0),
																																								ageOfMother(0),
																																								ageOfFather(0),
																																								tracked(0){
				// initialise two gametes for this individual
				Gamete gameteMaternal(p, rng);
				Gamete gametePaternal(p, rng);
																																												
				// fill the binary gene arrays
				geneticsBinary[0] = gameteMaternal.genesOfGamete;
				geneticsBinary[1] = gametePaternal.genesOfGamete;
																																												
				// fill the age-specific gene arrays
				ageSpecificGenesMaternal = gameteMaternal.ageSpecificGenesOfGamete;
				ageSpecificGenesPaternal = gametePaternal.ageSpecificGenesOfGamete;
				
				// calculates survivalProb, based on binary genes and fills averageSurvivalProbAgeGenes array
																																												// based on age-specific genes
				calcSurvivalProb(p);
				quality = p.initAgeSpecificGenes; // get initial quality

				// if the individual is female she should make gametes, otherwise the male should make stem cells.
				(isFemale) ? makeSeveralGametes(p, rng) : makeStemcells(p);
				isFemaleSex = (isFemale) ? 1 : 0;
}

Individual::Individual(Individual& mother,
																							Individual& father,
																							Randomizer& rng,
																							const Parameters& p) : age(0),
																																												  tracked(0),
																																												  ageOfMother(mother.age),
																																												  ageOfFather(father.age){
				/**Constructor to reproduce and create offspring . **/
																																																		
				// first, get a gamete from the mothers gamete list
				Gamete gameteMother = std::move(mother.gametes.back());
																																																		
				// remove this gamete from the mothers gamete list
				mother.gametes.pop_back();
																																																		
				// have the father generate a gamete
				Gamete gameteFather = father.makeGameteFromStemCell(rng);
				// TODO: have father duplicate the stemcell
				
				// make a new individual of these gametes
				geneticsBinary[0] = gameteMother.genesOfGamete;
				geneticsBinary[1] = gameteFather.genesOfGamete;
				ageSpecificGenesMaternal = gameteMother.ageSpecificGenesOfGamete;
				ageSpecificGenesPaternal = gameteFather.ageSpecificGenesOfGamete;
																																																		
				calcSurvivalProb(p); // to set the survival probability of the new individual
																																																		
				quality = averageSurvivalProbAgeGenes[age]; // get new individual its quality
																																																		
				// calculate effect of quality from both parents
				float parentalQuality = p.weightMaternalEffect * mother.quality + (1 - p.weightMaternalEffect) * father.quality;
																																																	
				if (p.addQuality)	survivalProb *= parentalQuality; // multiply survival prob with the quality of the parents
}

Gamete Individual::makeGamete(Randomizer& rng){
    /**Function to make a single gamete. Based on stochasticity to determine which genes the gamete receives. **/
				
    Gamete gamete; // initialise gamete
				
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // based on a random value, determine whether the genes of the gamete come from the mother or from the father
        gamete.genesOfGamete[i] = geneticsBinary[rng.bernoulli()][i];
    }
				
				for (size_t i = 0; i < gamete.ageSpecificGenesOfGamete.size(); ++i){
								gamete.ageSpecificGenesOfGamete[i] = (rng.bernoulli()) ? ageSpecificGenesMaternal[i] : ageSpecificGenesPaternal[i];
				}
				
    return gamete;
}

void Individual::makeSeveralGametes(const Parameters &p,
                             Randomizer& rng){
    /** Function to generate multiple gametes for a female.**/
				
    for (unsigned i = 0; i < p.numOfGametes; ++i){
        Gamete gamete = makeGamete(rng);
								gametes.push_back(gamete);
    }
}

bool Individual::dies(Randomizer& rng,
                      const Parameters& p){
    /**Function to calculate which individuals will die.**/
				
    bool dies = false;
    
				// multiply the survival probability of the individual at its age times
    // the effect of an extrinsic mortality risk on the survival probability
    const float survivalProbForAge = averageSurvivalProbAgeGenes[age] * survivalProb;
				// the survival prob based on both gene arrays
				
    float survivalProbIncExtrinsicRisk = survivalProbForAge * (1 - p.extrinsicMortRisk);
				// taking extrinsic mortality into account
				
    if (rng.bernoulli(survivalProbIncExtrinsicRisk)){ // bernoulli distribution with the bias of survival probability of the individual
        age += 1; // increment age if individual survives the mortality round
								quality = averageSurvivalProbAgeGenes[age]; // every time an individual ages, the parental quality is recalculated
        if (age == p.maximumAge) dies = true;
    } else { // indidvidual dies
        dies = true; // Individual will die
    }
    return dies;
}

void Individual::mutateGametes(const Parameters &p,
                               Randomizer &rng){
				
    for (size_t i = 0; i < gametes.size(); ++i){
								gametes[i].mutate(p, rng, false); // false refers to being a stem cell
    }
    
}

void Individual::makeStemcells(const Parameters& p){
				
    for (unsigned i = 0; i < p.numOfStemCells; ++i){
        Gamete gamete;
        gamete.genesOfGamete = geneticsBinary[0];
								gamete.ageSpecificGenesOfGamete = ageSpecificGenesMaternal;
        Gamete gamete2;
        gamete2.genesOfGamete = geneticsBinary[1];
								gamete2.ageSpecificGenesOfGamete = ageSpecificGenesPaternal;
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
				
    // first, get a random stem cell ...
				std::array<Gamete, 2> stemCell = stemCells[rng.drawRandomNumber(stemCells.size())];
				// ... then remove this stem cell from the list
				
    // initialise gamete
    Gamete gamete;
				// TODO: I use recombination in the generation of the gamete. Try out one array
				// draw one random number which is at least 40 bits long. Random long long int = 8 bytes
				// eg [0010101101010100010101001], n = 40. Pick every gene based on this
    
				for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // determine based on a bernoulli distribution which gene will be inherited
        gamete.genesOfGamete[i] = stemCell[rng.bernoulli()].genesOfGamete[i];
    }
				
				for (size_t i = 0; i < gamete.ageSpecificGenesOfGamete.size(); ++i){
								gamete.ageSpecificGenesOfGamete[i] = stemCell[rng.bernoulli()].ageSpecificGenesOfGamete[i];
				}
				
    return gamete;
}

void Individual::calcSurvivalProb(const Parameters& p){
				// calculate the survival probability based on the binary represented gene array
				// sum number of ones to calculate the survival probability
				int sumOfDamage1 = std::accumulate(geneticsBinary[0].begin(), geneticsBinary[0].end(), 0);
				int sumOfDamage2 = std::accumulate(geneticsBinary[1].begin(), geneticsBinary[1].end(), 0);
				
				
				// calculate survival probability based on number of ones
				survivalProb = exp(p.strengthOfSelection * (sumOfDamage1 + sumOfDamage2));
				
				// calculate the array with survival probabilities for the age-specific genes
				for (auto i = 0u; i < ageSpecificGenesMaternal.size(); ++i){
								float average = (ageSpecificGenesMaternal[i] + ageSpecificGenesPaternal[i]) * 0.5;
								averageSurvivalProbAgeGenes[i] = average;
				}
}
