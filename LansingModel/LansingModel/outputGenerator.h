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
        ofs << t << " "
        << i.age << " "
        << i.ageOfMother << " "
        << i.ageOfFather << " "
        << (i.averageAgeSpecificGenes[i.age] * i.survivalProb) << " "
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
								float ageSpecSurvProb; // male survival probability based on age-specific gene
								// if age-specific gene effect should not be taken into account, this will be set to 1.
								ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : male.averageAgeSpecificGenes[male.age];
								// calculates survival into the next year
								float s = ageSpecSurvProb * male.survivalProb * (1 - p.extrinsicMortRisk);
								// calculate expected age at death
        float expectedAgeAtDeath = male.age + (s / (1 - s));
        
        // write maternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfMother << " "
        << "F "
        << (ageSpecSurvProb * male.survivalProb) << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
        
        // write paternal information
        ofs << male.age << " "
        << expectedAgeAtDeath << " "
        << male.ageOfFather << " "
        << "M "
        << (ageSpecSurvProb * male.survivalProb) << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
    }

    
    for (Individual female : females){
								float ageSpecSurvProb; // female survival probability based on age-specific gene
								// if age-specific gene effect should not be taken into account, this will be set to 1.
								ageSpecSurvProb = (p.addQuality && !p.addAgeSpecific) ? 1 : female.averageAgeSpecificGenes[female.age];
								// calculates survival into the next year
        float s = ageSpecSurvProb * female.survivalProb * (1 - p.extrinsicMortRisk);
								// calculate expected age at death
        float expectedAgeAtDeath = female.age + (s / (1 - s));
        
        // write maternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfMother << " "
        << "F "
        << (ageSpecSurvProb * female.survivalProb) << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
        
        // write paternal information
        ofs << female.age << " "
        << expectedAgeAtDeath << " "
        << female.ageOfFather << " "
        << "M "
        << (ageSpecSurvProb * female.survivalProb) << " "
        << p.mutationProbStemcell << " "
        << p.mutationProb << std::endl;
    }
    ofs.close();
}

void createOutputTrackedIndividuals(const Parameters& p,
																																			const indVec& deadIndividuals){
				
				/**Function to create output for the tracked individuals. For every individual its age and the corresponding parental quality value
					is documented to the first file. To the second file information about the offspring of the tracked individuals is documented. **/
				
				// open file to write output for the
				std::ofstream ofs;
				ofs.open("outputWithParentalQuality.txt"); // the output file
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
												ofs << ind << " " // write as id of individual
												<< i << " " // write age
												<< deadIndividuals[ind].averageAgeSpecificGenes[i] << " " // write parental quality for this age class
												<< p.mutationProbAgeSpecificGenes << " "
												<< p.meanMutationBias << " "
												<< p.sdMutationalEffectSize << std::endl;
								}
								
								// write expected age at death of the offspring to a file
								for (size_t i = 0; i < deadIndividuals[ind].offspringOfIndividual.size(); ++i){ // loop through the number of offspring this individual has
												ofs2 << ind << " "; // use this index as identifier of the individual
												
												// calculate expected age at death for this offspring of the tracked individual
												// get the survival probability of the age-dependent genes
												float ageDependentSurvProb = deadIndividuals[ind].offspringOfIndividual[i].averageAgeSpecificGenes[deadIndividuals[ind].offspringOfIndividual[i].age];
												if (p.addQuality && !p.addAgeSpecific) ageDependentSurvProb = 1; // if
												// get the survival probability of the binary genes
												float binarySurvProb = deadIndividuals[ind].offspringOfIndividual[i].survivalProb;
												// mulitply the above-mentioned two to get total survival probability
												float totSurvProb = ageDependentSurvProb * binarySurvProb;
												// calculate s = yearly probability of survival to the next year
												float s = totSurvProb * (1 - p.extrinsicMortRisk);
												// calculate expected age at death
												float expectedAgeAtDeath = deadIndividuals[ind].offspringOfIndividual[i].age + (s / (1 - s));
												
												// if this flagged individual is male, the age of the father needs to be documented, if female > age of mother will be documented
												(deadIndividuals[ind].isFemaleSex) ? ofs2 << deadIndividuals[ind].offspringOfIndividual[i].ageOfMother : ofs2 << deadIndividuals[ind].offspringOfIndividual[i].ageOfFather;
												ofs2 << " ";
												(deadIndividuals[ind].isFemaleSex) ? ofs2 << "F " : ofs2 << "M ";
												ofs2 << (ageDependentSurvProb * binarySurvProb) << " " // write survival probability to file
												<< expectedAgeAtDeath << " "
												<< p.mutationProbAgeSpecificGenes << " "
												<< p.meanMutationBias << " "
												<< p.sdMutationalEffectSize << std::endl;
								}
				}
		
				ofs.close();
				ofs2.close();
}
