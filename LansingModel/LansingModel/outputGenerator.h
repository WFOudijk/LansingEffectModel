//
//  outputGenerator.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once

#include <fstream>
#include <cstdlib>
#include <vector>

using indVec = std::vector<Individual>;

void createOutput(const indVec& individualVec){
    /**Function to create output to the terminal of Individuals in a vector. **/
    std::cout << "Individual vector size = " << individualVec.size() << std::endl; 
    for (auto individual : individualVec){ // for every individual in the vector
        for (auto i = 0u; i < individual.genetics[0].genesOfGamete.size(); ++i) {
                std::cout << individual.genetics[0].genesOfGamete[i] << " ";
                std::cout << individual.genetics[1].genesOfGamete[i] << " ";
            }
        std::cout << std::endl;
        if (!individual.gametesOfIndividual.empty()){
            for (auto i : individual.gametesOfIndividual[2].genesOfGamete) std::cout << i << " ";
        }
        std::cout << std::endl;
        }
    }

void createOutputAgeDeath(const int t,
                          const Parameters& p,
                          const std::vector<int>& ageAtDeath){
    /**This function creates output to research the average death age. **/
    if (t == 0) { // if time is at zero. Empty file before addition
        std::ofstream ofs;
        ofs.open("outputDeathAge.csv"); // output file
        if (!ofs.is_open()){
            std::cerr << "Error. Unable to open output file.\n";
            exit(EXIT_FAILURE);
        }
        ofs.close();
    }
    std::ofstream ofs;
    ofs.open("outputDeathAge.csv", std::ios::app); // output file for age of death
    ofs << t << " "
    << p.mutationProb << " "
    << p.extrinsicMortRisk << " "
    << p.strengthOfSelection << " "
    << p.populationSize << " "
     << std::accumulate(ageAtDeath.begin(), ageAtDeath.end(), 0.0) /
        ageAtDeath.size() << std::endl; // look at average age of death over time
    ofs.close();
}

void createOutputDeclineInGameteQuality(const int t,
                                        const std::vector<Individual>& ageAtDeath){
    if (t == 0) { // if time is at zero. Empty file before addition
        std::ofstream ofs;
        ofs.open("outputDeclineGameteQuality.csv"); // output file
        if (!ofs.is_open()){
            std::cerr << "Error. Unable to open output file.\n";
            exit(EXIT_FAILURE);
        }
        ofs.close();
    }
    std::ofstream ofs;
    ofs.open("outputDeclineGameteQuality.csv", std::ios::app); // output file for age of death
    for (auto i : ageAtDeath){
        ofs << t << " "
        << i.age << " "
        << i.ageOfMother << " "
        << i.ageOfFather << std::endl; 
    }
    ofs.close();
}
