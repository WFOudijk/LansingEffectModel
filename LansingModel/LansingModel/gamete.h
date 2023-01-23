//
//  gamete.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include <array>

const int numOfGenes = 20; // the number of genes every individual contains
using arrayOfGenes = std::array<bool, numOfGenes>; // binary

struct Gamete{
    std::array<bool, numOfGenes> genesOfGamete; // an array with its genes
    // true (1) represents damage
    Gamete(){} // default constructor
    Gamete(const Parameters& p,
           Randomizer& rng){
        /** Constructor to initialise the gamete genetics.
        initDamageProportion represents the proportion of genes
        which need to be initialised to one, representing the initial damage. **/
        for (int i = 0; i < numOfGenes; ++i){ // to set some initial damage
            genesOfGamete[i] = rng.bernoulli(p.initDamageProportion);
        }
    }    
    void mutate(const Parameters& p, Randomizer& rng);
};

void Gamete::mutate(const Parameters &p,
                    Randomizer &rng){
    for (int i = 0u; i < genesOfGamete.size(); ++i){
        if (rng.bernoulli(p.mutationProb)){ // if mutation occurs:
            genesOfGamete[i] = 1; // the gene becomes damaged, meaning it becomes one.
        }
    }
   
}
