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

using arrayOfGenes = std::array<bool, numOfGenes>;
using vectorOfAgeSpecificGenes = std::vector<float>;

struct Individual {
    unsigned int age;
    unsigned int ageOfMother;
    unsigned int ageOfFather;
    float survivalProb; // based on the binary genes
    //float parentalQuality; // parental quality
    bool tracked; // the flagged individuals will be true
    bool isDead{false}; // track dead individuals
			
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
    void dies(Randomizer& rng, const Parameters& p);
    void mutateGametes(const Parameters& p, Randomizer& rng);
    void makeStemcells(const Parameters& p);
    void mutateStemCells(const Parameters& p, Randomizer& rng);
    Gamete makeGameteFromStemCell(const Parameters& p, Randomizer& rng);
    void calcSurvivalProb(const Parameters& p);
    void calcAverageParentalQuality();
    void calcAverageInvestmentGenes();
    unsigned int calcNumberOfOffspring(const Parameters& p, Randomizer& rng);
    void reproduce(const Parameters& p, Randomizer& rng, Individual& male);
    
};

Individual::Individual(const Parameters& p, // initializing constructor
                       Randomizer& rng,
                       bool isFemale) : age(0),
                                        ageOfMother(0), // should this not be >0 (for consistency)
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
   //parentalQuality = p.initAgeSpecificGenes;

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
                                                                                                                                                                                        
    //parentalQuality = averageAgeSpecificGenes[0]; // get new individual its quality
                                                                                                                                                                                        
    // calculate effect of quality from both parents
    //float effectQuality = p.weightMaternalEffect * mother.parentalQuality + (1.0 - p.weightMaternalEffect) * father.parentalQuality;
    float effectQuality = p.weightMaternalEffect *
                            mother.averageAgeSpecificGenes[mother.age] +
                            (1.0 - p.weightMaternalEffect) *
                            father.averageAgeSpecificGenes[father.age];

                                                  
    if (p.addQuality) survivalProb *= effectQuality; // multiply survival prob with the quality of the parents
                                                  
    // check if the investment genes affect the quality of the offspring
    if (p.addInvestmentAffectingOffspringQuality) {
    
        // get females investment in reproduction
        double investmentInReproduction = 1.0 - mother.averageInvestmentGenes[mother.age];
        
        // calculate the effect of the investment per offspring
        double adjustedInvestmentInReproduction = investmentInReproduction / p.numOfOffspringPerFemale;
        // TODO: is both investment effects are on, it should be divided by the number of offspring determined by the allocation effect
        
        // use the adjustedInvestmentInRepair to calculate the effect on the survival of the offspring
        survivalProb = survivalProb -
                ((p.baselineSurvival - adjustedInvestmentInReproduction) * p.scalingStrengthOfAllocationToSurvival);
        if (survivalProb < 0) survivalProb = 0; // to prevent the survival probability becoming negative
    }
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
    const std::bitset<41> y{rng.rui64()};
				
    for (size_t i = 0; i <= p.maximumAge; ++i){ // fill for every age class
        gamete.ageSpecificGenesOfGamete.push_back(ageSpecificGenes[y[i]][i]);
    }
    
    // get another random bitset, this needs to be the length of the maximum age for setting the age-specific genes
    const std::bitset<41> z{rng.rui64()};

    for (size_t i = 0; i <= p.maximumAge; ++i){ // fill for every age class
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
    /**Function to determine which individuals will die.**/
				
    // get age-specific survival probability
    float survivalProbAgeSpec = averageAgeSpecificGenes[age];
    // if quality is enabled but age-specific genetic effects are disabled, the age-specific genes do need to evolve (for quality determination)
    // but they should not be taken into account for determining survival probability of the individual
    if (p.addQuality && !p.addAgeSpecific) survivalProbAgeSpec = 1.0;
    
    // set investment in repair based on the individual's resource budget
    float investmentInRepair = 1.0 - p.weightInvestment * sqr(1.0 - averageInvestmentGenes[age]); // 1 - c3 * (1 - a)^2
    // if investment is disabled in the model. It should not play a part.
    if (!p.addInvestmentInRepair & !p.addInvestmentAffectingOffspringQuality) investmentInRepair = 1.0;
    //if (!p.addInvestmentInRepair) investmentInRepair = 1;

    // Caclulate the survival prob based on both gene arrays
    float totalSurvivalProb = survivalProbAgeSpec * survivalProb * investmentInRepair;
    // taking extrinsic mortality into account
    float adjustedSurvivalProb = totalSurvivalProb * p.survivalProbExtrinsicMort;
    
    if (rng.bernoulli(adjustedSurvivalProb)){ // bernoulli distribution with the bias of survival probability of the individual
        age += 1; // increment age if individual survives the mortality round
        //parentalQuality = averageAgeSpecificGenes[age]; // every time an individual ages, the parental quality is recalculated
        
        if (age == p.maximumAge) isDead = true;
    } else { // indidvidual dies
        isDead = true; // Individual will die.
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
    
    stemCells.reserve(p.numOfStemCells);
    
    stemCells.resize(p.numOfStemCells, genetics);

    
}

void Individual::mutateStemCells(const Parameters& p,
                                 Randomizer& rng){
    /**Function to mutate the male stem cells. Occurs every time step. **/
				
    for (auto& stemCell : stemCells){
        // true refers to them being stem cells
        stemCell[0].mutate(p, rng, true);
        stemCell[1].mutate(p, rng, true);
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
    // sum number of ones to calculate the survival probability
    auto sumOfDamage = std::accumulate(geneticsBinary[0].begin(), geneticsBinary[0].end(), 0) +
        std::accumulate(geneticsBinary[1].begin(), geneticsBinary[1].end(), 0);
				
    // calculate survival probability based on number of ones
    survivalProb = exp(p.strengthOfSelection * sumOfDamage);
}

void Individual::calcAverageParentalQuality(){
    /** Function to calculate the array with the average for the age-specific parental quality genes **/
    
    const auto& ageGenes1 = ageSpecificGenes[0];
    const auto& ageGenes2 = ageSpecificGenes[1];
    for (size_t i = 0; i < ageGenes1.size(); ++i){
        double averageGene = (ageGenes1[i] + ageGenes2[i]) * 0.5;
        averageAgeSpecificGenes.push_back(averageGene);
    }
}

void Individual::calcAverageInvestmentGenes(){
    /** Function to calculate the array with the average for the age-specific repair/reproduction distribution genes **/

    for (size_t i = 0; i < ageSpecificInvestmentInRepair[0].size(); ++i){
        averageInvestmentGenes.push_back((ageSpecificInvestmentInRepair[0][i] + ageSpecificInvestmentInRepair[1][i]) * 0.5);
    }
}

unsigned int Individual::calcNumberOfOffspring(const Parameters& p,
                                               Randomizer& rng){
    /**Function to calculate the number of offspring an individual wil get based on their investment in repair/reproduction distribution.
     The function uses stochastic rounding to determine the actual number of offspring. **/
    
    // determine investment in reproduction based on investment in repair genes.
    float investmentInReproduction = (1.0 - averageInvestmentGenes[age]);
    
    // calculate numOfOffspring based on sigmoidal/logistic function
    // 10 is to scale the numbers between 0 and 100 and to make the graph less steep.
    float numOfOffspring = p.scalingParameterForNumOfOffspring / (1.0 + exp(-p.scaleInvestmentValuesForCalc * (investmentInReproduction - p.pointOfHalfScalingParam))); // sigmoidal
    
    // stochastic rounding
    return stochasticRound(numOfOffspring, rng);

}

