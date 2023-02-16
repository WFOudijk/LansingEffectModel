//
//  gamete.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include <array>

const int numOfGenes = 20; // the number of genes every individual contains

struct Gamete{
    
    std::array<bool, numOfGenes> genesOfGamete; // array with binary genes
    // true (1) represents damage

    std::array<double, numOfGenes> ageSpecificGenesOfGamete; // array with age-specific genes
    // represented by survival probabilities
    
    int numOfMuts = 0;
    
    Gamete() : numOfMuts(0){} // default constructor
    
    Gamete(const Parameters& p,
           Randomizer& rng){
        /** Constructor to initialise the gamete genetics.
        initDamageProportion represents the proportion of genes
        which need to be initialised to one, representing the initial damage.
									the age-specific genes are initially filled with initAgeSpecificGenes. **/
        for (int i = 0; i < numOfGenes; ++i){ // to set some initial damage
            genesOfGamete[i] = rng.bernoulli(p.initDamageProportion);
        }
								ageSpecificGenesOfGamete.fill(p.initAgeSpecificGenes);
    }
    void mutate(const Parameters& p, Randomizer& rng, const bool isStemcell);
};

void Gamete::mutate(const Parameters &p,
                    Randomizer &rng,
                    const bool isStemcell){
				/**Function to mutate a gamete.  First a check occurs whether the gamete comes from a stem cell or from a female. Next, it checks for every gene
					if a mutation will occur. If so, for the age-specific genes the mutational effect is drawn from a normal distribution. If a mutation occurs in a binary gene
					this gene value will be set to one to indicate damage.
					**/
				for (size_t i = 0; i < genesOfGamete.size(); ++i){
								if (rng.bernoulli(p.mutationProbAgeSpecificGenes)){ // if mutation occurs:
												ageSpecificGenesOfGamete[i] += rng.drawMutationEffect(); // the gene is mutated based on distribution.
												if (ageSpecificGenesOfGamete[i] < 0) ageSpecificGenesOfGamete[i] = 0; // negative survival probability makes no sense
												if (ageSpecificGenesOfGamete[i] > 1) ageSpecificGenesOfGamete[i] = 1; // survival probability bigger than 1 makes no sense
								}
								if (isStemcell) { // if the stemcell will mutate
												if (rng.bernoulli(p.mutationProbStemcell)) genesOfGamete[i] = 1; // make this gene damaged
        } else { // else a gamete will mutate
												// check seperately for the gene array with binary genes if a mutation occurs.
												if (rng.bernoulli(p.mutationProb)) genesOfGamete[i] = 1; // make this gene damaged
        }
    }
}
