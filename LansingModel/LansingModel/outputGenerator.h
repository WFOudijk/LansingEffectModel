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

void createOutputDeclineInGameteQuality(const int t,
                                        const Parameters& p,
                                        const std::vector<Individual>& deadIndividualsVec){
    
    /**Function to write information about the death individuals every p.outputTime step**/
    
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
    ofs.open("outputDeclineGameteQuality.txt", std::ios::app); // output file
    for (auto i : deadIndividualsVec){
        // if age-specific gene effect should not be taken into account, this will be set to 1.
        float ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : i.averageAgeSpecificGenes[i.age];
        // set investment in repair to 1
        float investmentInRepair = 1.0;
        // if investment is enabled in the model. It should play a part.
        if (p.addInvestmentInRepair) investmentInRepair = 1 - p.weightInvestment * sqr(1 - i.averageInvestmentGenes[i.age]); // 1 - c3 * (1 - a)^2
        
        char sex = i.isFemaleSex ? 'F' : 'M';
        ofs << t << " "
        << i.age << " " // is the age at death of individual
        << sex << " "
        << i.ageOfMother << " "
        << i.ageOfFather << " "
        << (ageSpecSurvProb * i.survivalProb * investmentInRepair) << " "
        << i.averageInvestmentGenes[i.age] << "\n";
    }
    ofs.close();
}

void createOutputAgeSpecificGenes(const Parameters& p,
                                    const indVec& deadIndividuals){
				
    /**Function to write the gene values for every age class to a file. Function is called after time simulation.
     So the final gene values are documented. **/
				
    // open file to write output for the
    std::ofstream ofs;
    ofs.open("outputWithAgeSpecificGenes.txt");
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }

    // write parental quality for every age class to file
    for (size_t ind = 0; ind < deadIndividuals.size(); ++ind){
        // write to first file
        for (unsigned i = 0; i < p.maximumAge; ++i){
            ofs << ind << "_" << p.mutationProb << "_" << p.mutationProbStemcell << "_"
            << p.meanMutationBias << "_" << p.sdMutationalEffectSize << "_"
            << p.mutationProbAgeSpecificGenes << "_" << p.weightInvestment << " " // write as id of individual
            << i << " " // write age
            << deadIndividuals[ind].averageInvestmentGenes[i] << " " // write parental investment in repair for this age class
            << deadIndividuals[ind].averageAgeSpecificGenes[i] << "\n"; // write parental quality for this age class
        }
		
    }
    ofs.close();
}

void outputForSimulatedLifespan(const indVec& deadIndividuals){
    /**Function to write relevant output of the offspring to determine life expectancy
     for the cross-sectional analysis. **/

    // open file to write output for the
    std::ofstream ofs;
    ofs.open("outputLifeExp.txt");
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }
    
    // loop through the dead offspring
    for (size_t ind = 0; ind < deadIndividuals.size(); ++ind){
        ofs << ind << " " // use index as ID
        << deadIndividuals[ind].age << " " // age at death of individual
        << deadIndividuals[ind].ageOfMother << " " // write maternal age at birth
        << deadIndividuals[ind].ageOfFather << "\n"; // write paternal age at birth
    }

    ofs.close();
}

void outputOffspringLifespanLongitudinal(indVec deadTrackedIndividuals){
    /**Function to write output of the longitudinal offspring lifespan simulation. **/
    
    // open file to write output for the longitudinal offspring lifespan simulation
    std::ofstream ofs;
    ofs.open("outputLifeExpLong.txt");
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }
    
    for (size_t ind = 0; ind < deadTrackedIndividuals.size(); ++ind){
        for (Individual offspring : deadTrackedIndividuals[ind].offspring){
            ofs << ind << " " // use this as ID of the parent
            << offspring.age << " " // the age at death of this individual
            << offspring.ageOfMother << " " // write maternal age to file
            << offspring.ageOfFather << "\n"; // write paternal age to file
        }
    }
    ofs.close();
}
