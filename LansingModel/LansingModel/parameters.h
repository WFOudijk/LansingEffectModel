//
//  parameters.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//
#pragma once
#include <fstream>

struct Parameters {
    // initialise the parameters
    Parameters() : populationSize(100),
                   initDamageProportion(0.5),
                   numOfOffspringPerFemale(2),
                   mutationProb(0.001), // based on parameter simulations
                   extrinsicMortRisk(0.05), // maximum added number of years to live will be 19
                   outputTime(10),
                   tEnd(10000), // 10.000
                   strengthOfSelection(-0.05),
                   maximumAge(40),
                   mutationProbStemcell(0.001),
                   meanMutationBias(-0.02),
                   sdMutationalEffectSize(0.01),
                   initAgeSpecificGenes(0.99),
																   mutationProbAgeSpecificGenes(0.0001), // with 0.005 all three models are statistically significant. When only ageSpecific is true
                   numOfIndividualsToFollow(100),
																   addBinary(true),
																   addAgeSpecific(true),
																   addQuality(true){ // based on parameter simulations
                       numOfGametes = maximumAge * numOfOffspringPerFemale;
                       numOfStemCells = numOfGametes * 2;
                   }
    
    unsigned int populationSize, numOfOffspringPerFemale, maximumAge, numOfIndividualsToFollow, numOfGametes,numOfStemCells; // total population size
    //float initDamageProportion; // the proportion of initial damage in the genome
    //unsigned int numOfOffspringPerFemale; // number of offspring a female should produce
    //float mutationProb; // probability a mutation will occur
    //float extrinsicMortRisk; // the extrinsic mortality risk, equal for every adult
    int outputTime; // when to output info
    int tEnd; // end of simulation
    //float strengthOfSelection; // this coefficient determines the strength of the effect of damage
    //unsigned int maximumAge; // maximum age an individual can get to
    float mutationProbStemcell, meanMutationBias, sdMutationalEffectSize, initAgeSpecificGenes, initDamageProportion, mutationProb, extrinsicMortRisk, strengthOfSelection, mutationProbAgeSpecificGenes; // mutation probability of stemcell to mutate
    //float meanMutationBias;
    //float sdMutationalEffectSize;
    //float initAgeSpecificGenes;
				//float mutationProbAgeSpecificGenes; // mutation probability for the age specific genes
				//unsigned int numOfIndividualsToFollow, numOfGametes, numOfStemCells, populationSize, numOfOffspringPerFemale, maximumAge; // number of individuals to follow longitudinal
				//unsigned int numOfGametes; // the derived number of gametes a female should have
				//unsigned int numOfStemCells; // number of stem cells a male should create
				bool addBinary, addAgeSpecific, addQuality;
				//bool addAgeSpecific;
				//bool addQuality;
    
    void readParameters(const std::string& parameterFile);
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    float& parameter,
                    std::ifstream& ifs);
				void checkParam(const std::string parID,
																				const std::string focal_parametername,
																				unsigned int& parameter,
																				std::ifstream& ifs);
};

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            float& parameter,
                            std::ifstream& ifs) {
    // set parameter from file to parameter in object parameter
    if (parID == focal_parametername) {
        ifs >> parameter;
        std::clog << "Parameter " << parID << " is set to " << parameter << std::endl;
    }
}

void Parameters::checkParam(const std::string parID,
																												const std::string focal_parametername,
																												unsigned int& parameter,
																												std::ifstream& ifs) {
				// set parameter from file to parameter in object parameter
				if (parID == focal_parametername) {
								ifs >> parameter;
								std::clog << "Parameter " << parID << " is set to " << parameter << std::endl;
				}
}

void Parameters::readParameters(const std::string& parameterFile){
    /**This function receives a parameter file and reads this. Next, the parameters in the file are set to the correct parameters
     in the parameters object using the checkParam function.**/
    std::ifstream ifs(parameterFile.c_str());
    if(!ifs.is_open()){
        std::cerr << "Error. Unable to read the following parameter file: "
                    << parameterFile << std::endl;
        exit(EXIT_FAILURE);
    }
    std::clog << "Reading parameters from file: " << parameterFile << std::endl;
    for(;;){
        std::string parID;
        ifs >> parID; // get row in file
        if(ifs.good()) { // setting of the parameters
            checkParam(parID, "mutationProb", mutationProb, ifs);
												checkParam(parID, "mutationProbStemcell", mutationProbStemcell, ifs);
            checkParam(parID, "mutationProbAgeSpecificGenes", mutationProbAgeSpecificGenes, ifs);
												checkParam(parID, "numOfIndividualsToFollow", numOfIndividualsToFollow, ifs);
												checkParam(parID, "populationSize", populationSize, ifs);
        }
        else break;
    }
}

