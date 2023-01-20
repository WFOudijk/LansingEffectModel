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
    Parameters() : populationSize(1000), // 500, 1500, 2000 > parameter exploration
                   initDamageProportion(0.5),
                   numOfOffspringPerFemale(2),
                   mutationProb(0.01), // 0.01, 0.15, steps: 0.02. > parameter exploration
                   extrinsicMortRisk(0.05),
                   outputTime(5),
                   tEnd(10000), // 10.000
                   maximumAge(40){}
    
    int populationSize; // total population size
    double halfPopulation; // to determine number of males and females
    double initDamageProportion;
    int numOfOffspringPerFemale; // number of offspring a female should produce
    double mutationProb; // probability a mutation will occur
    double extrinsicMortRisk; // the extrinsic mortality risk, equal for every adult
    int outputTime; // when to output info
    int tEnd; // end of simulation
    int maximumAge; 
    
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
            //checkParam(parID, "extrinsicMortRisk", extrinsicMortRisk, ifs);
        }
        else break;
    }
}

