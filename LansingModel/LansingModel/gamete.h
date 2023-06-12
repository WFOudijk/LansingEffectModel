//
//  gamete.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include "parameters.h"
#include "randomizer.h"
#include "utils.h"
#include <array>

const int numOfGenes = 20; // the number of genes every individual contains
using arrayOfGenes = std::array<bool, numOfGenes>;
using vectorOfAgeSpecificGenes = std::vector<float>;

struct Gamete{

    std::array<bool, numOfGenes> genesOfGamete; // array with binary genes
    // true (1) represents damage

    std::vector<float> ageSpecificGenesOfGamete; // array with age-specific genes
    // represented by floats between 0 and 1

    std::vector<float> ageSpecificInvestmentInRepair; // array with age-specific investment in repair values

    Gamete(){} // default constructor

    Gamete(const Parameters& p,
           Randomizer& rng){
				
        /** Constructor to initialise the gamete genetics.
         initDamageProportion represents the proportion of genes
         which need to be initialised to one, representing the initial damage.
         the age-specific genes are initially filled with initAgeSpecificGenes. **/
				
        for (int i = 0; i < numOfGenes; ++i){ // to set some initial damage
            genesOfGamete[i] = rng.bernoulli(p.initDamageProportion);
        }

        // sets age-specific array to length maximum age and filled with intial value
        ageSpecificGenesOfGamete.resize((p.maximumAge + 1), p.initAgeSpecificGenes);

        // set age-specific invstment in repair to maximum age and fill with initial value
        ageSpecificInvestmentInRepair.resize((p.maximumAge + 1), p.initInvestmentInRepair);
    }

    void mutate(const Parameters& p, Randomizer& rng, const bool isStemcell);
};

void Gamete::mutate(const Parameters &p,
                    Randomizer &rng,
                    const bool isStemcell){

    /**Function to mutate a gamete.  First is the age-specific gene array mutated.
     This is done by using a Poisson distribution to draw how many times a mutation will occur. Next, the genes are randomly sampled to determine which will mutate.
     Next, the binary gene array is mutated.
     This is done by first checking whether the gamete comes from a male stem cell or from a female gamete. Next, another Poisson distribution is used
     to determine how many times a mutation will occur. Finally, the genes are randomly sampled to determine which will mutate.
     **/

    // if both are off, the genes don't need to mutate. Saves speed in model.
    if(p.addQuality || p.addAgeSpecific){
        // mutation of age-specific genes
        int numOfEvents = rng.drawNumOfMuts(); // draws how many mutations will occur
        // mutate
        for (int i = 0; i < numOfEvents; ++i){
            // determine which gene will mutate
            int geneToMutate = rng.drawRandomNumber(ageSpecificGenesOfGamete.size());
            // draw the effect from the mutation based on a normal distribution
            ageSpecificGenesOfGamete[geneToMutate] += rng.drawMutationEffect();
            // check if the gene value is not < 0 or > 1. If so, gene value is clipped
            clip01(ageSpecificGenesOfGamete[geneToMutate]);
        }
    }
    
    if (p.addInvestmentInRepair) { // if resrouce distribution is on, these genes will mutate
        // mutation of age-specific genes for investment in repair/ reproduction
        const double expectedNumMut{ageSpecificInvestmentInRepair.size() * p.mutationProbInvestmentGenes};
        const unsigned numMut{rng.rpois(expectedNumMut)};
        for (unsigned i = 0; i < numMut; ++i){
            int geneToMutate = rng.drawRandomNumber(ageSpecificInvestmentInRepair.size());
            ageSpecificInvestmentInRepair[geneToMutate] += rng.drawMutationEffectInvestment();
            clip01(ageSpecificInvestmentInRepair[geneToMutate]);
        }
    }
    
    // mutation of binary genes
    // get mutation probability depending on if it is a stem cell or gamete mutating
    double mutationProb = isStemcell ? p.mutationProbStemcell : p.mutationProb;
    
    // get number of mutation based on poisson distribution
    unsigned numMut{rng.rpois(genesOfGamete.size() * mutationProb)};
    
    // mutate
    for (size_t i = 0; i < numMut; ++i){
        auto test = rng.rn(genesOfGamete.size());
        genesOfGamete[test] = 1;
    }
}
