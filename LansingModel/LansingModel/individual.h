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
#include <cassert>

struct Individual {
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    float survivalProb; // based on the binary genes
    bool tracked; // to keep track of individuals for longitudinal analysis
    bool isDead{false}; // to track dead individuals
			
    // array with genes - binary    
    std::array<arrayOfGenes, 2> geneticsBinary;

    // array with age specific genes - survival probabilities/ parental care qualities
    std::array<vectorOfAgeSpecificGenes, 2> ageSpecificGenes;
							    
    // averages of the maternal and paternal survival probabilities
    vectorOfAgeSpecificGenes averageAgeSpecificGenes;
				
    // array with age-specific investment in repair vs reproduction genes
    std::array<vectorOfAgeSpecificGenes, 2> ageSpecificInvestmentInRepair;
    
    // averages of the age-specific investment genes
    vectorOfAgeSpecificGenes averageInvestmentGenes;
    
    // vector with gametes if individual is female
    std::vector<Gamete> gametes;
    // vector with stem cells if individual is male
    std::vector<std::array<Gamete, 2> > stemCells;
    
    // vector with offspring of individual for the longitudinal analysis
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
    void dies(Randomizer& rng, const Parameters& p);
    void mutateGametes(const Parameters& p, Randomizer& rng);
    void makeStemcells(const Parameters& p);
    void mutateStemCells(const Parameters& p, Randomizer& rng);
    Gamete makeGameteFromStemCell(const Parameters& p, Randomizer& rng);
    void calcSurvivalProb(const Parameters& p);
    void calcAverageParentalQuality();
    void calcAverageInvestmentGenes();
    void reproduce(const Parameters& p, Randomizer& rng, Individual& male);
    
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
																																												
    // fill the age-specific survival gene arrays
    ageSpecificGenes[0] = gameteMaternal.ageSpecificGenesOfGamete;
    ageSpecificGenes[1] = gametePaternal.ageSpecificGenesOfGamete;
																																																																																		
    // fill the age-specific investment in repair arrays
    ageSpecificInvestmentInRepair[0] = gameteMaternal.ageSpecificInvestmentInRepair;
    ageSpecificInvestmentInRepair[1] = gametePaternal.ageSpecificInvestmentInRepair;

    // calculates survivalProb, based on binary genes
    calcSurvivalProb(p);
    // calculates and fills averageAgeSpecificGenes array
    calcAverageParentalQuality();
    // calculates and fills averageInvestmentGenes
    calcAverageInvestmentGenes();
    // simulate parental quality effect
    if (p.addQuality) survivalProb *= averageAgeSpecificGenes[0];
   

    // if the individual is female she should make gametes, otherwise the male should make stem cells.
    (isFemale) ? makeSeveralGametes(p, rng) : makeStemcells(p);
    isFemaleSex = isFemale;
}

Individual::Individual(Individual& mother,
                       Individual& father,
                       Randomizer& rng,
                       const Parameters& p) : age(0),
                                              ageOfMother(mother.age),
                                              ageOfFather(father.age),
                                              tracked(0){
    /**Constructor to create a new offspring because of reproduction. **/
																																																		
    // first, get a gamete from the mothers gamete list
    Gamete gameteMother = std::move(mother.gametes.back());
                                                                                                                                                                                        
    // remove this gamete from the mothers gamete list
    mother.gametes.pop_back();
                                                                                                                                                                                        
    // have the father generate a gamete from a stem cell
    Gamete gameteFather = father.makeGameteFromStemCell(p, rng);

    // make a new individual of these gametes
    geneticsBinary[0] = gameteMother.genesOfGamete;
    geneticsBinary[1] = gameteFather.genesOfGamete;
    ageSpecificGenes[0] = gameteMother.ageSpecificGenesOfGamete;
    ageSpecificGenes[1] = gameteFather.ageSpecificGenesOfGamete;
    ageSpecificInvestmentInRepair[0] = gameteMother.ageSpecificInvestmentInRepair;
    ageSpecificInvestmentInRepair[1] = gameteFather.ageSpecificInvestmentInRepair;
    
    calcSurvivalProb(p); // to set the survival probability of the new individual
    calcAverageParentalQuality(); // averages the age-specific survival gene arrays
    calcAverageInvestmentGenes(); // average the age-specific investment gene arrays;
    
    if (p.addQuality) { // if quality is enabled in the model ..
        // .. calculate effect of quality from both parents
        float effectQuality = p.weightMaternalEffect *
                                mother.averageAgeSpecificGenes[mother.age] +
                                (1.0 - p.weightMaternalEffect) *
                                father.averageAgeSpecificGenes[father.age];
        
        // multiply survival prob of the new individual with the quality effect of the parents
        survivalProb *= effectQuality;
    }
    
    // check if the investment genes affect the quality of the offspring
    if (p.addInvestmentInRepair) {
    
        // get females investment in reproduction
        double investmentInReproduction = 1.0 - mother.averageInvestmentGenes[mother.age];
        
        // calculate the effect of the investment per offspring
        double adjustedInvestmentInReproduction = investmentInReproduction / p.numOfOffspringPerFemale;
        
        // use the adjustedInvestmentInRepair to calculate the effect on the survival of the offspring
        survivalProb *= logistic(p.steepnessAllocationToReproduce, adjustedInvestmentInReproduction, p.scalingStrengthOfAllocationToReproduce);
    }
}

Gamete Individual::makeGamete(Randomizer& rng,
                              const Parameters& p){
    /**Function to make a single gamete. Based on recombination to determine which genes the gamete receives. **/
				
    Gamete gamete; // initialise gamete
				
    // get a bit of numOfGenes long to use the bits as random for setting the gamete genes
    const std::bitset<numOfGenes> x{rng.rui32()};
				
    for (int i = 0; i < numOfGenes; ++i){ // loop through every gene
        // based on the random bit sequence, determine if gene is paternal or maternal
        gamete.genesOfGamete[i] = geneticsBinary[x[i]][i];
    }
				
    // get another random bitset, this needs to be the length of the maximum age for setting the age-specific genes
    const std::bitset<41> y{rng.rui64()};
				
    for (size_t i = 0; i <= p.maximumAge; ++i){ // fill for every age class
        // again, using the bit sequence determine if gene is paternal or maternal
        gamete.ageSpecificGenesOfGamete.push_back(ageSpecificGenes[y[i]][i]);
    }
    
    // get another random bitset, this needs to be the length of the maximum age for setting the age-specific genes
    const std::bitset<41> z{rng.rui64()};

    for (size_t i = 0; i <= p.maximumAge; ++i){ // fill for every age class
        // again, using the bit sequence, determine if gene is paternal or maternal
        gamete.ageSpecificInvestmentInRepair.push_back(ageSpecificInvestmentInRepair[z[i]][i]);
    }
				
    return gamete;
}

void Individual::makeSeveralGametes(const Parameters &p,
                                    Randomizer& rng){
    /** Function to generate multiple gametes for a female.**/
				
    for (unsigned i = 0; i < p.numOfGametes; ++i){
        Gamete gamete = makeGamete(rng, p);
        gametes.emplace_back(gamete);
    }
}

void Individual::dies(Randomizer& rng,
                      const Parameters& p){
    /**Function to determine if the individual will die.**/
				
    // set age-specific survival probability to 1 to disable it
    float survivalProbAgeSpec = 1.0;
    // if age-specific genetic effect is enabled, the age-specific genes determine (a.o.) survival
    if (p.addAgeSpecific) survivalProbAgeSpec = averageAgeSpecificGenes[age];
    
    // set investment in repair to 1 to disable it
    float investmentInRepair = 1.0;
    // if investment is enabled in the model. It should be recalculated based on the genes
    if (p.addInvestmentInRepair) investmentInRepair = 1.0 - p.weightInvestment * sqr(1.0 - averageInvestmentGenes[age]); // 1 - c3 * (1 - a)^2

    // Caclulate the survival prob based on both gene arrays and the binary survival effect
    float totalSurvivalProb = survivalProbAgeSpec * survivalProb * investmentInRepair;
    // taking extrinsic mortality into account
    float adjustedSurvivalProb = totalSurvivalProb * p.survivalProbExtrinsicMort;
    
    // draw number from bernoulli distribution based on survival probability of the individual
    if (rng.bernoulli(adjustedSurvivalProb)){ // individual survives
        age += 1; // increment age
        if (age == p.maximumAge) isDead = true; // if the individual has reached maximum age, they will die
    } else { // indidvidual dies
        isDead = true;
    }
}

void Individual::mutateGametes(const Parameters &p,
                               Randomizer &rng){
    /**Function to mutate the gametes of a female. **/
				
    for (auto& gamete : gametes){
        gamete.mutate(p, rng, false); // false refers to being a stem cell
    }
    
}

void Individual::makeStemcells(const Parameters& p){
    /**Function to make stem cells for a male. These stem cells are represented by gametes and have an
     identical genome as the individual.  **/
				
    Gamete gamete;
    gamete.genesOfGamete = geneticsBinary[0];
    gamete.ageSpecificGenesOfGamete = ageSpecificGenes[0];
    gamete.ageSpecificInvestmentInRepair = ageSpecificInvestmentInRepair[0];
    Gamete gamete2;
    gamete2.genesOfGamete = geneticsBinary[1];
    gamete2.ageSpecificGenesOfGamete = ageSpecificGenes[1];
    gamete2.ageSpecificInvestmentInRepair = ageSpecificInvestmentInRepair[1];
    
    std::array<Gamete, 2> genetics = {gamete, gamete2};
    
    stemCells.reserve(p.numOfStemCells);
    
    stemCells.resize(p.numOfStemCells, genetics);
    
}

void Individual::mutateStemCells(const Parameters& p,
                                 Randomizer& rng){
    /**Function to mutate the male stem cells. **/
				
    for (auto& stemCell : stemCells){
        // true refers to them being stem cells
        stemCell[0].mutate(p, rng, true);
        stemCell[1].mutate(p, rng, true);
    }
}

Gamete Individual::makeGameteFromStemCell(const Parameters& p,
                                          Randomizer& rng){
    /**Function to make a gamete using a male stem cell. This will occur when a male is selected to reproduce.
     Recombination occurs here. */
				
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
    const std::bitset<41> y{rng.rui64()};
				
    for (size_t i = 0; i <= p.maximumAge; ++i){
        // for every age-specific gene of new gamete, determine which gene is inherited
        gamete.ageSpecificGenesOfGamete.push_back(stemCell[y[i]].ageSpecificGenesOfGamete[i]);
    }
    
    // get anther template random bitset, this one the length of maximumage to determine age-specific gene distribution
    const std::bitset<41> z{rng.rui64()};

    for (size_t i = 0; i <= p.maximumAge; ++i){
        // for every age-specific gene of new gamete, determine which gene is inherited
        gamete.ageSpecificInvestmentInRepair.push_back(stemCell[z[i]].ageSpecificInvestmentInRepair[i]);
    }
				
    return gamete;
}

void Individual::calcSurvivalProb(const Parameters& p){
    /**Function to calculate the survival probability of the individual.
     Based on the number of damaged genes in the binary array. */
    
    // calculate the survival probability based on the binary represented gene array
    // sum the number of ones (so the number of damaged genes) to calculate the survival probability
    auto sumOfDamage = std::accumulate(geneticsBinary[0].begin(), geneticsBinary[0].end(), 0) +
        std::accumulate(geneticsBinary[1].begin(), geneticsBinary[1].end(), 0);
				
    // calculate survival probability based on number of ones
    survivalProb = exp(p.strengthOfSelection * sumOfDamage);
}

void Individual::calcAverageParentalQuality(){
    /** Function to calculate the average of the set of age-specific survival genes. **/
    
    for (size_t i = 0; i < ageSpecificGenes[0].size(); ++i){
        // fill the corresponding array with the averages
        averageAgeSpecificGenes.push_back((ageSpecificGenes[0][i] + ageSpecificGenes[1][i]) * 0.5);
    }
}

void Individual::calcAverageInvestmentGenes(){
    /** Function to calculate the average of the set of age-specific repair/reproduction genes **/

    for (size_t i = 0; i < ageSpecificInvestmentInRepair[0].size(); ++i){
        averageInvestmentGenes.push_back((ageSpecificInvestmentInRepair[0][i] + ageSpecificInvestmentInRepair[1][i]) * 0.5);
    }
}

