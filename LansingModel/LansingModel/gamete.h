//
//  gamete.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include <array>

const int numOfGenes = 20; // the number of genes every individual contains
//using arrayOfGenes = std::array<bool, numOfGenes>; // binary

struct Gamete{
    
    std::array<double, numOfGenes> genesOfGamete; // an array with its genes
    
    int numOfMuts = 0;
    // true (1) represents damage
    
    Gamete() : numOfMuts(0){} // default constructor
    
    Gamete(const Parameters& p,
           Randomizer& rng){
        /** Constructor to initialise the gamete genetics.
        initDamageProportion represents the proportion of genes
        which need to be initialised to one, representing the initial damage. **/
        for (int i = 0; i < numOfGenes; ++i){ // to set some initial damage
            genesOfGamete[i] = rng.bernoulli(p.initDamageProportion);
        }
    }
//    // move constructor
//    Gamete(Gamete&& other) : genesOfGamete(std::move(other.genesOfGamete)),
//                             numOfMuts(std::move(other.numOfMuts)) {}
    void mutate(const Parameters& p, Randomizer& rng, const bool isStemcell);
};

void Gamete::mutate(const Parameters &p,
                    Randomizer &rng,
                    const bool isStemcell){
    if (isStemcell) { // if the stemcell will mutate
        for (size_t i = 0; i < genesOfGamete.size(); ++i){
            if (rng.bernoulli(p.mutationProbStemcell)){ // if mutation occurs:
                genesOfGamete[i] += rng.drawMutationEffect(); // the gene is mutated based on distribution.
                if (genesOfGamete[i] < 0) genesOfGamete[i] = 0; // negative survival probability makes no sense
                if (genesOfGamete[i] > 1) genesOfGamete[i] = 1; // survival probability bigger than 1 makes no sense
            }
        }
    } else { // else a gamete will mutate
        for (size_t i = 0; i < genesOfGamete.size(); ++i){
            if (rng.bernoulli(p.mutationProb)){ // if mutation occurs:
                numOfMuts += 1;
                genesOfGamete[i] += rng.drawMutationEffect(); // the gene is mutated based on distribution.
                if (genesOfGamete[i] < 0) genesOfGamete[i] = 0; // negative survival probability makes no sense
                if (genesOfGamete[i] > 1) genesOfGamete[i] = 1; // survival probability bigger than 1 makes no sense
            }
        }
    }
}
