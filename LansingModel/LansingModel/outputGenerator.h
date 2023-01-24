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
#include "calculateAverage.h"

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

//void createOutputForPlot(const indVec& males,
//                         const indVec& females){
//    /**This function creates the output for a scatter plot where age is at the x-axis and average survival probability at the y-axis. **/
//    std::ofstream ofs;
//    ofs.open("output.csv"); // output file
//    if (!ofs.is_open()){
//        std::cerr << "Error. Unable to open output file.\n";
//        exit(EXIT_FAILURE);
//    }
//    // calculate average survival prob per age of the females
//    auto female_avg = calcAverageAcrossAgeClasses(females);
//    // calculate average survival prob per age of the males
//    auto male_avg = calcAverageAcrossAgeClasses(males);
//
//    for (auto i = 0u; i < female_avg.size(); ++i) { // loop through every age
//        double pop_avg = (female_avg[i] + male_avg[i]) * 0.5; // calculate population average
//        ofs << pop_avg << std::endl; // write this average to file
//    }
//    ofs.close();
//}
//
//void createOuputForGGPlot(const indVec& males,
//                          const indVec& females,
//                          const int t,
//                          const Parameters& p){
//    /**Function to create the output for a GGPlot. It determines the average survival probability for both the females and males,
//     then it determines the population average based on this. This is written to a file, including the current time and the age. **/
//    if (t == 0) { // if time is at zero. Empty file before addition
//        std::ofstream ofs;
//        ofs.open("outputFacetWrap.csv"); // output file
//        if (!ofs.is_open()){
//            std::cerr << "Error. Unable to open output file.\n";
//            exit(EXIT_FAILURE);
//        }
//        ofs.close();
//    }
//
//    std::ofstream ofs;
//    ofs.open("outputFacetWrap.csv", std::ios::app); // output file
//    if (!ofs.is_open()){
//        std::cerr << "Error. Unable to open output file.\n";
//        exit(EXIT_FAILURE);
//    }
//    auto femaleAverage = calcAverageAcrossAgeClasses(females); // calculate average survival prob per age of the females
//    auto maleAverage = calcAverageAcrossAgeClasses(males); // calculate average survival prob per age of the males
//
//    for (auto i = 0u; i < femaleAverage.size(); ++i) { // loop through every age
//        double popAverage = (femaleAverage[i] + maleAverage[i]) * 0.5; // calculate population average
//        ofs << t << " "
//            << p.meanMutationBias << " "
//            << p.sdMutationalEffectSize << " "
//            << p.extrinsicMortRisk << " "
//            << i << " "
//            << popAverage << std::endl; // write current time, age and average to file
//    }
//    ofs.close();
//}
//
//void createOutputLifeExpectancy(const indVec& males,
//                                const indVec& females,
//                                const Parameters& p){
//    /**This function can be used to create output to look at the individual life expectancy. **/
//    std::vector<double> malesLE = calcLifeExpectancyPerIndividual(males);
//    std::vector<double> femalesLE = calcLifeExpectancyPerIndividual(females);
//    std::ofstream ofs;
//    ofs.open("outputLE.csv"); // output file
//    if (!ofs.is_open()){
//        std::cerr << "Error. Unable to open output file.\n";
//        exit(EXIT_FAILURE);
//    }
//    for (auto i = 0u; i < malesLE.size(); ++i){
//        ofs << p.meanMutationBias << " "
//        << p.sdMutationalEffectSize << " "
//        << p.extrinsicMortRisk << " "
//        << malesLE[i] << " " << femalesLE[i] << std::endl;
//    }
//    ofs.close();
//}

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
