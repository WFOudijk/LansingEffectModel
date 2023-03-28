//
//  individual.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include "gamete.h"

#include <vector>
#include <bitset>

using arrayOfGenes = std::array<bool, numOfGenes>;
using vectorOfAgeSpecificGenes = std::vector<float>;

struct Individual {
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    float survivalProb; // based on the binary genes
    float parentalQuality; // parental quality
    bool tracked; // the flagged individuals will be true
			
    // array with genes - binary
    std::array<arrayOfGenes, 2> geneticsBinary;

    // array with age specific genes - survival probabilities
    std::array<vectorOfAgeSpecificGenes, 2> ageSpecificGenes;
							    
    // averaging the maternal and paternal survival probabilities
    vectorOfAgeSpecificGenes averageAgeSpecificGenes;
				
    // array with age-specific investment in repair vs reproduction
    std::array<vectorOfAgeSpecificGenes, 2> ageSpecificInvestmentInRepair;
    
    // averaging the age-specific invetsment genes
    vectorOfAgeSpecificGenes averageInvestmentGenes;
    
    std::vector<Gamete> gametes;
    std::vector<std::array<Gamete, 2> > stemCells;
    
    std::vector<Individual> offspring;
    bool isFemaleSex;
    
    Individual(const Parameters& p, // initializing constructor
               Randomizer& rng,
               bool isFemale);

    Individual(Individual& mother,
               Individual& father,
               Randomizer& rng,
               const Parameters& p);
				
    Gamete makeGamete(Randomizer& rng, const Parameters& p);
    void makeSeveralGametes(const Parameters &p, Randomizer& rng);
    bool dies(Randomizer& rng, const Parameters& p);
    void mutateGametes(const Parameters& p, Randomizer& rng);
    void makeStemcells(const Parameters& p);
    void mutateStemCells(const Parameters& p, Randomizer& rng);
    Gamete makeGameteFromStemCell(const Parameters& p, Randomizer& rng);
    void calcSurvivalProb(const Parameters& p);
    void calcAverageParentalQuality();
    void calcAverageInvestmentGenes();
    unsigned int calcNumberOfOffspring(const Parameters& p, Randomizer& rng);
    
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
    ageSpecificGenes[0] = gameteMaternal.ageSpecificGenesOfGamete;
    ageSpecificGenes[1] = gametePaternal.ageSpecificGenesOfGamete;
																																																																																		
    // fill the age-specific investment in repair array
    ageSpecificInvestmentInRepair[0] = gameteMaternal.ageSpecificInvestmentInRepair;
    ageSpecificInvestmentInRepair[1] = gametePaternal.ageSpecificInvestmentInRepair;

    // calculates survivalProb, based on binary genes and
    calcSurvivalProb(p);
    // fills averageAgeSpecificGenes array based on age-specific genes
    calcAverageParentalQuality();
    // fills average averageInvestmentGenes
    calcAverageInvestmentGenes();
    // get initial quality
   parentalQuality = p.initAgeSpecificGenes;

    // if the individual is female she should make gametes, otherwise the male should make stem cells.
    (isFemale) ? makeSeveralGametes(p, rng) : makeStemcells(p);
    isFemaleSex = (isFemale) ? 1 : 0;
}

Individual::Individual(Individual& mother,
                       Individual& father,
                       Randomizer& rng,
                       const Parameters& p) : age(0),
                                              ageOfMother(mother.age),
                                              ageOfFather(father.age),
                                              tracked(0){
    /**Constructor to reproduce and create offspring . **/
																																																		
    // first, get a gamete from the mothers gamete list
    Gamete gameteMother = std::move(mother.gametes.back());
                                                                                                                                                                                        
    // remove this gamete from the mothers gamete list
    mother.gametes.pop_back();
                                                                                                                                                                                        
    // have the father generate a gamete
    Gamete gameteFather = father.makeGameteFromStemCell(p, rng);

    // make a new individual of these gametes
    geneticsBinary[0] = gameteMother.genesOfGamete;
    geneticsBinary[1] = gameteFather.genesOfGamete;
    ageSpecificGenes[0] = gameteMother.ageSpecificGenesOfGamete;
    ageSpecificGenes[1] = gameteFather.ageSpecificGenesOfGamete;
    ageSpecificInvestmentInRepair[0] = gameteMother.ageSpecificInvestmentInRepair;
    ageSpecificInvestmentInRepair[1] = gameteFather.ageSpecificInvestmentInRepair;
    
    calcSurvivalProb(p); // to set the survival probability of the new individual
    calcAverageParentalQuality(); // averages the age-specific gene arrays
    calcAverageInvestmentGenes(); // average the age-specific investment gene arrays;
                                                                                                                                                                                        
    parentalQuality = averageAgeSpecificGenes[0]; // get new individual its quality
                                                                                                                                                                                        
    // calculate effect of quality from both parents
    float effectQuality = p.weightMaternalEffect * mother.parentalQuality + (1 - p.weightMaternalEffect) * father.parentalQuality;
                                                                                                                                                                                    
    if (p.addQuality) survivalProb *= effectQuality; // multiply survival prob with the quality of the parents
}

Gamete Individual::makeGamete(Randomizer& rng,
                              const Parameters& p){
    /**Function to make a single gamete. Based on stochasticity to determine which genes the gamete receives. **/
				
    Gamete gamete; // initialise gamete
				
    // get a bit of numOfGenes long to use the bits as random for setting the gamete genes
    const std::bitset<numOfGenes> x{rng.rui32()};
				
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // based on the random bit sequence, determine whether the genes of the gamete come from the mother or from the father
        gamete.genesOfGamete[i] = geneticsBinary[x[i]][i];
    }
				
    // get another random bitset, this needs to be the length of the maximum age for setting the age-specific genes
    const std::bitset<40> y{rng.rui64()};
				
    for (size_t i = 0; i < p.maximumAge; ++i){ // fill for every age class
        gamete.ageSpecificGenesOfGamete.push_back(ageSpecificGenes[y[i]][i]);
    }
    
    // get another random bitset, this needs to be the length of the maximum age for setting the age-specific genes
    const std::bitset<40> z{rng.rui64()};

    for (size_t i = 0; i < p.maximumAge; ++i){ // fill for every age class
        gamete.ageSpecificInvestmentInRepair.push_back(ageSpecificInvestmentInRepair[z[i]][i]);
    }
				
    return gamete;
}

void Individual::makeSeveralGametes(const Parameters &p,
                                    Randomizer& rng){
    /** Function to generate multiple gametes for a female.**/
				
    for (unsigned i = 0; i < p.numOfGametes; ++i){
        Gamete gamete = makeGamete(rng, p);
        gametes.push_back(gamete);
    }
}

bool Individual::dies(Randomizer& rng,
                      const Parameters& p){
    /**Function to determine which individuals will die.**/
				
    bool dies = false;
    
    // get age-specific survival probability
    float survivalProbAgeSpec = averageAgeSpecificGenes[age];
    // if quality is on but age-specific genetic effects is off, the age-specific genes do need to evolve (for quality determination)
    // but they should not be taken into account for determining survival probability of the individual
    if (p.addQuality && !p.addAgeSpecific) survivalProbAgeSpec = 1;
    
    // set investment in repair based on their resource budget
    float investmentInRepair = 1.0 - p.weightInvestment * sqr(1.0 - averageInvestmentGenes[age]); // 1 - c3 * (1 - a)^2
    
    // if investment is off in the model. It should not play a part. 
    if (!p.addInvestmentInRepair) investmentInRepair = 1;
    // the survival prob based on both gene arrays
    float survivalProbTot = survivalProbAgeSpec * survivalProb * investmentInRepair;
    // taking extrinsic mortality into account
    float survivalProbIncExtrinsicRisk = survivalProbTot * (1 - p.extrinsicMortRisk);
    
    if (rng.bernoulli(survivalProbIncExtrinsicRisk)){ // bernoulli distribution with the bias of survival probability of the individual
        age += 1; // increment age if individual survives the mortality round
        parentalQuality = averageAgeSpecificGenes[age]; // every time an individual ages, the parental quality is recalculated
        if (age == p.maximumAge) dies = true;
    } else { // indidvidual dies
        dies = true; // Individual will die
    }
    return dies;
}

void Individual::mutateGametes(const Parameters &p,
                               Randomizer &rng){
    /**Function to mutate the gametes of a female. **/
				
    for (size_t i = 0; i < gametes.size(); ++i){
        gametes[i].mutate(p, rng, false); // false refers to being a stem cell
    }
    
}

void Individual::makeStemcells(const Parameters& p){
    /**Function to make stem cells for a male. These stem cells are represented by gametes and have the same genome as the
     male whose stem cells these are.  **/
				
    Gamete gamete;
    gamete.genesOfGamete = geneticsBinary[0];
    gamete.ageSpecificGenesOfGamete = ageSpecificGenes[0];
    gamete.ageSpecificInvestmentInRepair = ageSpecificInvestmentInRepair[0];
    Gamete gamete2;
    gamete2.genesOfGamete = geneticsBinary[1];
    gamete2.ageSpecificGenesOfGamete = ageSpecificGenes[1];
    gamete2.ageSpecificInvestmentInRepair = ageSpecificInvestmentInRepair[1];
    std::array<Gamete, 2> genetics = {gamete, gamete2};
				
    for (unsigned i = 0; i < p.numOfStemCells; ++i){
        stemCells.push_back(genetics);
    }
}

void Individual::mutateStemCells(const Parameters& p,
                                 Randomizer& rng){
    /**Function to mutate the male stem cells. Occurs every time step. **/
				
    for (size_t i = 0; i < stemCells.size(); ++i){
        // true refers to them being stem cells 
        stemCells[i][0].mutate(p, rng, true);
        stemCells[i][1].mutate(p, rng, true);
    }
}

Gamete Individual::makeGameteFromStemCell(const Parameters& p,
                                          Randomizer& rng){
    /**Function to make a gamete using a male stem cell. This will occur when a male is selected to reproduce. */
				
    // first, get a random stem cell
    std::array<Gamete, 2> stemCell = stemCells[rng.drawRandomNumber(stemCells.size())];
				
    // initialise gamete
    Gamete gamete;
    
    // draw a random numOfGenes long bitset to use as a random template for distributing genes
    const std::bitset<numOfGenes> x{rng.rui32()};
				
    for (int i = 0; i < numOfGenes; ++i){ // loop through every binary gene
        // determine based on the bitset which binary gene will be inherited
        gamete.genesOfGamete[i] = stemCell[x[i]].genesOfGamete[i];
    }
				
    // get anther template random bitset, this one the length of maximumage to determine age-specific gene distribution
    const std::bitset<40> y{rng.rui64()};
				
    for (size_t i = 0; i < p.maximumAge; ++i){
        // for every age-specific gene of new gamete, determine which gene is inherited
        gamete.ageSpecificGenesOfGamete.push_back(stemCell[y[i]].ageSpecificGenesOfGamete[i]);
    }
    
    // get anther template random bitset, this one the length of maximumage to determine age-specific gene distribution
    const std::bitset<40> z{rng.rui64()};

    for (size_t i = 0; i < p.maximumAge; ++i){
        // for every age-specific gene of new gamete, determine which gene is inherited
        gamete.ageSpecificInvestmentInRepair.push_back(stemCell[z[i]].ageSpecificInvestmentInRepair[i]);
    }
				
    return gamete;
}

void Individual::calcSurvivalProb(const Parameters& p){
    /**Function to calculate the survival probability of the individual. Based on the number of damaged genes in the binary array. */
    
    // calculate the survival probability based on the binary represented gene array
    // sum number of ones to calculate the survival probability
    int sumOfDamage1 = std::accumulate(geneticsBinary[0].begin(), geneticsBinary[0].end(), 0);
    int sumOfDamage2 = std::accumulate(geneticsBinary[1].begin(), geneticsBinary[1].end(), 0);
				
    // calculate survival probability based on number of ones
    survivalProb = exp(p.strengthOfSelection * (sumOfDamage1 + sumOfDamage2));
}

void Individual::calcAverageParentalQuality(){
    /** Function to calculate the array with the average for the age-specific parental quality genes **/
    
    for (size_t i = 0; i < ageSpecificGenes[0].size(); ++i){
        averageAgeSpecificGenes.emplace_back((ageSpecificGenes[0][i] + ageSpecificGenes[1][i]) * 0.5);
    }
}

void Individual::calcAverageInvestmentGenes(){
    /** Function to calculate the array with the average for the age-specific repair/reproduction distribution genes **/

    for (size_t i = 0; i < ageSpecificInvestmentInRepair[0].size(); ++i){
        averageInvestmentGenes.emplace_back((ageSpecificInvestmentInRepair[0][i] + ageSpecificInvestmentInRepair[1][i]) * 0.5);
    }
}

unsigned int Individual::calcNumberOfOffspring(const Parameters& p,
                                               Randomizer& rng){
    /**Function to calculate the number of offspring an individual wil get based on their investment in repair/reproduction distribution.
     The function uses stochastic rounding to determine the actual number of offspring. **/
    
    // determine investment in reproduction based on investment in repair genes.
    float investmentInReproduction = (1 - averageInvestmentGenes[age]);
    
    // calculate numOfOffspring based on sigmoidal/logistic function
    // 10 is to scale the numbers between 0 and 100 and to make the graph less steep.
    //float numOfOffspring = p.maxOffspring / (1 + std::exp(-10 * (investmentInReproduction - p.pointOfHalfMaxOffspring))); // sigmoidal
    
    float numOfOffspring = p.maxOffspring * investmentInReproduction / (p.pointOfHalfMaxOffspring + investmentInReproduction); // quickest?
    
    //float numOfOffspring = p.maxOffspring * pow(investmentInReproduction, (p.pointOfHalfMaxOffspring + investmentInReproduction)); max * x^(0.5+x)
    
    // stochastic rounding
    float decimal = numOfOffspring - trunc(numOfOffspring); // get decimal
    // draw bernoulli based on decimal, if true: round up, if false: round down
    unsigned int numOfOffspringRounded = (rng.bernoulli(decimal)) ? (numOfOffspring + (1 - decimal)) : trunc(numOfOffspring);
    // return rounded number
    return numOfOffspringRounded;
}
