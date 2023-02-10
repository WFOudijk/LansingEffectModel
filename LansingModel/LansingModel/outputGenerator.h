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
        ofs.open("outputDeathAge.txt"); // output file
        if (!ofs.is_open()){
            std::cerr << "Error. Unable to open output file.\n";
            exit(EXIT_FAILURE);
        }
        ofs.close();
    }
    std::ofstream ofs;
    ofs.open("outputDeathAge.txt", std::ios::app); // output file for age of death
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
                                        const Parameters& p,
                                        const std::vector<Individual>& deadIndividualsVec){
    if (t == 0) { // if time is at zero. Empty file before addition
        std::ofstream ofs;
        ofs.open("outputDeclineGameteQuality.txt"); // output file
        if (!ofs.is_open()){
            std::cerr << "Error. Unable to open output file.\n";
            exit(EXIT_FAILURE);
        }
        ofs.close();
    }
    std::ofstream ofs;
    ofs.open("outputDeclineGameteQuality.txt", std::ios::app); // output file for age of death
    for (auto i : deadIndividualsVec){
        ofs << t << " "
        << i.age << " "
        << i.ageOfMother << " "
        << i.ageOfFather << " "
        << i.survivalProb << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
    }
    ofs.close();
}

void createOutputLifeExpectancy(const Parameters& p,
                                const indVec& males,
                                const indVec& females){
    std::ofstream ofs;
    ofs.open("outputLifeExpectancy.txt"); // the output file
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }
    for (Individual male : males){
        double s = male.survivalProb * (1 - p.extrinsicMortRisk);
        double expectedAgeAtDeath = male.age + (s / (1 - s));
        
        // write maternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfMother << " "
        << "F "
        << male.survivalProb << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
        
        // write paternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfFather << " "
        << "M "
        << male.survivalProb << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
    }

    
    for (Individual female : females){
        double s = female.survivalProb * (1 - p.extrinsicMortRisk);
        double expectedAgeAtDeath = female.age + (s / (1 - s));
        
        // write maternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfMother << " "
        << "F "
        << female.survivalProb << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
        
        // write paternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfFather << " "
        << "M "
        << female.survivalProb << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
    }
    ofs.close();
}

void createOutputTrackedIndividuals(const Parameters& p,
                                    const indVec& trackedDeadIndividuals){
    /**
     Function to write the output of the tracked individuals. These are longitudinally followed by keeping track of the offspring they
     get over their lifetime.
     */
    std::ofstream ofs;
    ofs.open("outputLETrackedIndividuals.txt"); // the output file
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }
    for (size_t ind = 0; ind < trackedDeadIndividuals.size(); ++ind){ // loop through the flagged individuals
        for (size_t i = 0; i < trackedDeadIndividuals[ind].offspringOfIndividual.size(); ++i){ // loop through the number of offspring this individual has
            ofs << ind << " "; // use this index as identifier of the individual
            
            // calculate expected age at death for this offspring of the tracked individual
            double s = trackedDeadIndividuals[ind].offspringOfIndividual[i].survivalProb * (1 - p.extrinsicMortRisk);
            double expectedAgeAtDeath = trackedDeadIndividuals[ind].offspringOfIndividual[i].age + (s / (1 - s));
            
            // if this flagged individual is male, the age of the father needs to be documented, if female > age of mother will be documented
            (trackedDeadIndividuals[ind].sex == 'M') ? ofs << trackedDeadIndividuals[ind].offspringOfIndividual[i].ageOfFather : ofs << trackedDeadIndividuals[ind].offspringOfIndividual[i].ageOfMother;
            ofs << " "
            << trackedDeadIndividuals[ind].offspringOfIndividual[i].survivalProb << " " // write survival probability to file
            << expectedAgeAtDeath << std::endl; // write expected age at death to file
        }
    }
    ofs.close();
}

