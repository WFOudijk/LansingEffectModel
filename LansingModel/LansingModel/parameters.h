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
                   mutationProb(0.0045), // based on parameter simulations
                   extrinsicMortRisk(0.05), // maximum added number of years to live will be 19
                   outputTime(10),
                   tEnd(10000), // 10.000
                   strengthOfSelection(-0.05),
                   maximumAge(40),
                   mutationProbStemcell(0.0045),
                   numOfIndividualsToFollow(500){ // based on parameter simulations
                       numOfGametes = maximumAge * numOfOffspringPerFemale;
                       numOfStemCells = numOfGametes * 2;
                   }
    
    unsigned int populationSize; // total population size
    double initDamageProportion; // the proportion of initial damage in the genome
    unsigned int numOfOffspringPerFemale; // number of offspring a female should produce
    double mutationProb; // probability a mutation will occur
    double extrinsicMortRisk; // the extrinsic mortality risk, equal for every adult
    int outputTime; // when to output info
    int tEnd; // end of simulation
    double strengthOfSelection; // this coefficient determines the strength of the effect of damage
    unsigned int maximumAge; // maximum age an individual can get to
    unsigned int numOfGametes; // the derived number of gametes a female should have
    unsigned int numOfStemCells; // number of stem cells a male should create
    double mutationProbStemcell; // mutation probability of stemcell to mutate
    unsigned int numOfIndividualsToFollow; // number of individuals to follow longitudinal
    
    void readParameters(const std::string& parameterFile);
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    double& parameter,
                    std::ifstream& ifs);
};

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            double& parameter,
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
            checkParam(parID, "extrinsicMortRisk", extrinsicMortRisk, ifs);
            checkParam(parID, "strengthOfSelection", strengthOfSelection, ifs);
            checkParam(parID, "mutationProbStemcell", mutationProbStemcell, ifs);

        }
        else break;
    }
}

