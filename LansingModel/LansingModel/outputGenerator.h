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
        // if age-specific gene effect should not be taken into account, this will be set to 1.
        float ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : i.averageAgeSpecificGenes[i.age];
        // set investment in repair based on their resource budget
        float investmentInRepair = 1 - p.weightInvestment * sqr(1 - i.averageInvestmentGenes[i.age]); // 1 - c3 * (1 - a)^2
        // if investment is off in the model. It should not play a part.
        if (!p.addInvestmentInRepair & !p.addInvestmentAffectingOffspringQuality) investmentInRepair = 1;
        
        char sex = i.isFemaleSex ? 'F' : 'M';
        ofs << t << " "
        << i.age << " "
        << sex << " "
        << i.ageOfMother << " "
        << i.ageOfFather << " "
        << (ageSpecSurvProb * i.survivalProb * investmentInRepair) << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << " "
        << p.mutationProbInvestmentGenes << " "
        << p.sdMutationalEffectInvestmentInRepair << " "
        << i.averageInvestmentGenes[i.age] << std::endl;
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
        // if age-specific gene effect should not be taken into account, this will be set to 1.
        float ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : male.averageAgeSpecificGenes[male.age];
        // set investment in repair based on their resource budget
        float investmentInRepair = 1 - p.weightInvestment * sqr(1 - male.averageInvestmentGenes[male.age]); // 1 - c3 * (1 - a)^2
        // if investment is off in the model. It should not play a part.
        if (!p.addInvestmentInRepair & !p.addInvestmentAffectingOffspringQuality) investmentInRepair = 1;
        
        
        // calculates survival into the next year
        float s = ageSpecSurvProb * male.survivalProb * investmentInRepair * (1 - p.extrinsicMortRisk);
        // calculate expected age at death
        float expectedAgeAtDeath = male.age + (s / (1 - s));
        
        // write maternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfMother << " "
        << "F "
        << (ageSpecSurvProb * male.survivalProb * investmentInRepair) << " "
        << p.mutationProb << " "
        << p.mutationProbStemcell << " "
        << p.meanMutationBias << " "
        << p.sdMutationalEffectSize << " "
        << p.mutationProbAgeSpecificGenes << std::endl;
        
        // write paternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfFather << " "
        << "M "
        << (ageSpecSurvProb * male.survivalProb * investmentInRepair) << " "
        << p.mutationProb << " "
        << p.mutationProbStemcell << " "
        << p.meanMutationBias << " "
        << p.sdMutationalEffectSize << " "
        << p.mutationProbAgeSpecificGenes << std::endl;
    }

    
    for (Individual female : females){
        float ageSpecSurvProb; // female survival probability based on age-specific gene
        // if age-specific gene effect should not be taken into account, this will be set to 1.
        ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : female.averageAgeSpecificGenes[female.age];
        // set investment in repair based on their resource budget
        float investmentInRepair = 1 - p.weightInvestment * sqr(1 - female.averageInvestmentGenes[female.age]); // 1 - c3 * (1 - a)^2
        // if investment is off in the model. It should not play a part.
        if (!p.addInvestmentInRepair & !p.addInvestmentAffectingOffspringQuality) investmentInRepair = 1;
        
        // calculates survival into the next year
        float s = ageSpecSurvProb * female.survivalProb * investmentInRepair * (1 - p.extrinsicMortRisk);
        // calculate expected age at death
        float expectedAgeAtDeath = female.age + (s / (1 - s));
        
        // write maternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfMother << " "
        << "F "
        << (ageSpecSurvProb * female.survivalProb * investmentInRepair) << " "
        << p.mutationProb << " "
        << p.mutationProbStemcell << " "
        << p.meanMutationBias << " "
        << p.sdMutationalEffectSize << " "
        << p.mutationProbAgeSpecificGenes << std::endl;
        
        // write paternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfFather << " "
        << "M "
        << (ageSpecSurvProb * female.survivalProb * investmentInRepair) << " "
        << p.mutationProb << " "
        << p.mutationProbStemcell << " "
        << p.meanMutationBias << " "
        << p.sdMutationalEffectSize << " "
        << p.mutationProbAgeSpecificGenes << std::endl;
    }
    ofs.close();
}

void createOutputTrackedIndividuals(const Parameters& p,
                                    const indVec& deadIndividuals){
				
    /**Function to create output for the tracked individuals. For every individual its age and the corresponding parental quality value
     is documented to the first file. To the second file information about the offspring of the tracked individuals is documented. **/
				
    // open file to write output for the
    std::ofstream ofs;
    //ofs.open("outputWithInvestmentDistribution.txt"); // the output file
    ofs.open("outputWithInvestment.txt");
    if (!ofs.is_open()){
        std::cerr << "Error. Unable to open output file.\n";
        exit(EXIT_FAILURE);
    }
				
				
    // open file to write output of the tracked individuals
    std::ofstream ofs2;
    ofs2.open("outputLETrackedIndividuals.txt"); // the output file
    if (!ofs2.is_open()){
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
            ///<< deadIndividuals[ind].averageAgeSpecificGenes[i] << " " // write parental quality in repair for this age class
            << p.mutationProb << " "
            << p.mutationProbStemcell << " "    
            << p.meanMutationBias << " "
            << p.sdMutationalEffectSize << " "
            << p.mutationProbAgeSpecificGenes << " "
            << p.mutationProbInvestmentGenes << " "
            << p.sdMutationalEffectInvestmentInRepair << std::endl;
        }
								
        // write expected age at death of the offspring to a file
        for (size_t i = 0; i < deadIndividuals[ind].offspring.size(); ++i){ // loop through the number of offspring this individual has
            ofs2 << ind << "_" << p.mutationProb << "_" << p.mutationProbStemcell << "_"
            << p.meanMutationBias << "_" << p.sdMutationalEffectSize << "_"
            << p.mutationProbAgeSpecificGenes << "_" << p.weightInvestment << " "; // use this index as identifier of the individual
            
            // calculate expected age at death for this offspring of the tracked individual
            // get the survival probability of the age-dependent genes
            float ageDependentSurvProb = deadIndividuals[ind].offspring[i].averageAgeSpecificGenes[deadIndividuals[ind].offspring[i].age];
            if (p.addQuality && !p.addAgeSpecific) ageDependentSurvProb = 1; // if
            
            // get the survival probability of the binary genes
            float binarySurvProb = deadIndividuals[ind].offspring[i].survivalProb;
            
            // set investment in repair based on their resource budget
            float investmentInRepair = 1 - p.weightInvestment * sqr(1 - deadIndividuals[ind].offspring[i].averageInvestmentGenes[deadIndividuals[ind].offspring[i].age]); // 1 - c3 * (1 - a)^2
            // if investment is off in the model. It should not play a part.
            if (!p.addInvestmentInRepair & !p.addInvestmentAffectingOffspringQuality) investmentInRepair = 1;
            
            // mulitply the above-mentioned three to get total survival probability
            float totSurvProb = ageDependentSurvProb * binarySurvProb * investmentInRepair;
            // calculate s = yearly probability of survival to the next year
            float s = totSurvProb * (1 - p.extrinsicMortRisk);
            // calculate expected age at death
            float expectedAgeAtDeath = deadIndividuals[ind].offspring[i].age + (s / (1 - s));
            // if this flagged individual is male, the age of the father needs to be documented, if female > age of mother will be documented
            (deadIndividuals[ind].isFemaleSex) ? ofs2 << deadIndividuals[ind].offspring[i].ageOfMother : ofs2 << deadIndividuals[ind].offspring[i].ageOfFather;
            ofs2 << " ";
            (deadIndividuals[ind].isFemaleSex) ? ofs2 << "F " : ofs2 << "M ";
            ofs2 << (ageDependentSurvProb * binarySurvProb * investmentInRepair) << " " // write survival probability to file
            << expectedAgeAtDeath << " "
            << p.mutationProb << " "
            << p.mutationProbStemcell << " "
            << p.meanMutationBias << " "
            << p.sdMutationalEffectSize << " "
            << p.mutationProbAgeSpecificGenes << " "
            << p.mutationProbInvestmentGenes << " "
            << p.sdMutationalEffectInvestmentInRepair << std::endl;
        }
    }
		
    ofs.close();
    ofs2.close();
}
