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
    Parameters() : populationSize(100), // N in manuscript
                   initDamageProportion(0.1),
                   numOfOffspringPerFemale(1), // o in manuscript
                   mutationProb(0.0024), // mu_b in manuscript
                   extrinsicMortRisk(0.0),
                   outputTime(10),
                   tEnd(1000),
                   strengthOfSelection(-0.05), // s in manuscript
                   maximumAge(40), // M in manuscript
                   mutationProbStemcell(0.0024), // mu_b in manuscript
                   meanMutationBias(-0.02), // b_s in manuscript
                   sdMutationalEffectSize(0.02), // sigma in manuscript
                   initAgeSpecificGenes(0.90),
                   mutationProbAgeSpecificGenes(0.002), // mu_a in manuscript
                   weightMaternalEffect(0.5),
                   initInvestmentInRepair(0.5),
                   numOfStemCells(30), // n_sc in manuscript
                   meanMutationBiasInvestmentInRepair(0), // b_r in manuscript
                   sdMutationalEffectInvestmentInRepair(0.02), // sigma in manuscript
                   mutationProbInvestmentGenes(0.002), // mu_a in manuscript
                   weightInvestment(0.3), // c in manuscript
                   scalingStrengthOfAllocationToReproduce(1), // d in manuscript
                   numOfOffspringForOffspringLifespanSim(10),
                   steepnessAllocationToReproduce(3), // a in manuscript
                   addBinary(true),
                   addAgeSpecific(false),
                   addQuality(false),
                   addInvestmentInRepair(false){}
    
    unsigned int populationSize; // total population size
    float initDamageProportion; // the proportion of initial damage in the binary genes
    unsigned int numOfOffspringPerFemale; // number of offspring a female should produce
    float mutationProb; // probability a mutation will occur
    float extrinsicMortRisk; // the extrinsic mortality risk, equal for every adult
    int outputTime; // when to output info
    unsigned int tEnd; // end of simulation
    float strengthOfSelection; // this coefficient determines the strength of the effect of damage
    unsigned int maximumAge; // maximum age an individual can get to
    float mutationProbStemcell; // mutation probability of stemcell to mutate
    float meanMutationBias; // the mean mutation bias of the age-specific genes
    float sdMutationalEffectSize; // the sd of the mutation bias, i.e., the effect size
    float initAgeSpecificGenes; // the initial value for the age-specific genes
    float mutationProbAgeSpecificGenes; // mutation probability for the age specific genes
    unsigned int numOfGametes; // the derived number of gametes a female should have
    float weightMaternalEffect; // to determine how much the maternal quality affects the offspring in ratio to paternal effect.
    float initInvestmentInRepair; // initial investement in repair vs reproduction
    unsigned int numOfStemCells; // number of stem cells a male should create
    float meanMutationBiasInvestmentInRepair; // mean for normal distribution to draw mutation effect on investment in repair genes
    float sdMutationalEffectInvestmentInRepair; // sd for normal distribution to draw mutation effect on investment in repair genes
    float mutationProbInvestmentGenes; // mutation rate of age-specific investment in repair genes
    float weightInvestment; // to weigh the investment in repair genes
    float scalingStrengthOfAllocationToReproduce; // to determine the strength of the allocative effect on survival
    int numOfOffspringForOffspringLifespanSim; // the number of offspring per female to track to determine offspring lifespan in cross-sectional analysis
    float survivalProbExtrinsicMort; // the survival probability based on extrinsic mortality
    float steepnessAllocationToReproduce; // to determine the steepness in allocation to survival effect
    bool addBinary; // add binary genes to model
    bool addAgeSpecific; // adds age-specific genes to model
    bool addQuality; // adds quality effect to model
    bool addInvestmentInRepair; // adds investment in repair to model
    
    void readParameters(const std::string& parameterFile);
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    float& parameter,
                    std::ifstream& ifs);
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    unsigned int& parameter,
                    std::ifstream& ifs);
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    bool& parameter,
                    std::ifstream& ifs);
    void setAdditionalParams();
};

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            float& parameter,
                            std::ifstream& ifs) {
    // set parameter from file to parameter in object float parameter
    if (parID == focal_parametername) {
        ifs >> parameter;
        std::clog << "Parameter " << parID << " is set to " << parameter << std::endl;
    }
}

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            unsigned int& parameter,
                            std::ifstream& ifs) {
    // set parameter from file to parameter in object int parameter
    if (parID == focal_parametername) {
        ifs >> parameter;
        std::clog << "Parameter " << parID << " is set to " << parameter << std::endl;
    }
}

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            bool& parameter,
                            std::ifstream& ifs) {
    // set parameter from file to parameter in object bool parameter
    if (parID == focal_parametername) {
        int ifs2;
        ifs >> ifs2;
        bool ifs3 = (ifs2 == 0 ? false : true);
        parameter = ifs3;
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
            checkParam(parID, "populationSize", populationSize, ifs);
            checkParam(parID, "meanMutationBias", meanMutationBias, ifs);
            checkParam(parID, "sdMutationalEffectSize", sdMutationalEffectSize, ifs);
            checkParam(parID, "mutationProbInvestmentGenes", mutationProbInvestmentGenes, ifs);
            checkParam(parID, "sdMutationalEffectInvestmentInRepair", sdMutationalEffectInvestmentInRepair, ifs);
            checkParam(parID, "weightInvestment", weightInvestment, ifs);
            checkParam(parID, "tEnd", tEnd, ifs);
            checkParam(parID, "addBinary", addBinary, ifs);
            checkParam(parID, "addAgeSpecific", addAgeSpecific, ifs);
            checkParam(parID, "addQuality", addQuality, ifs);
            checkParam(parID, "addInvestmentInRepair", addInvestmentInRepair, ifs);
        }
        else break;
    }
}

void Parameters::setAdditionalParams(){
    /**Function to set some additional parameters. This is done after the parameter file is read. **/
    
    // every female can get, at most, maximumAge * numOfOffspringPerFemale number of offspring for the time loop
    // for the offspring lifespan simulation every female gets numOfOffspringForOffspringLifespanSim number of
    // offspring. Total number of gametes necessary is thus the following:
    numOfGametes = maximumAge * numOfOffspringPerFemale + numOfOffspringForOffspringLifespanSim;
    
    // set the survival probability based on the extrinsic mortality risk
    survivalProbExtrinsicMort = 1.0 - extrinsicMortRisk;
}

