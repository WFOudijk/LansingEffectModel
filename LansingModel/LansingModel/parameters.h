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
    Parameters() : populationSize(1000),
                   initDamageProportion(0.5),
                   numOfOffspringPerFemale(3),
                   mutationProb(0.0024),
                   extrinsicMortRisk(0.05), // maximum added number of years to live will be 19
                   outputTime(10),
                   tEnd(200000), // 10.000
                   strengthOfSelection(-0.05),
                   maximumAge(40),
                   mutationProbStemcell(0.0024),
                   meanMutationBias(-0.022),
                   sdMutationalEffectSize(0.024),
                   initAgeSpecificGenes(0.99),
                   mutationProbAgeSpecificGenes(0.003), // 0.003 based on param simulations
                   numOfIndividualsToFollow(500),
                   weightMaternalEffect(0.5),
                   initInvestmentInRepair(0.5),
                   numOfStemCells(30),
                   meanMutationBiasInvestmentInRepair(0),
                   sdMutationalEffectInvestmentInRepair(0.02),
                   mutationProbInvestmentGenes(0.004),
                   weightInvestment(0.3),
                   scalingParameterForNumOfOffspring(4),
                   pointOfHalfScalingParam(0.5),
                   baselineSurvival(0.5),
                   scalingStrengthOfAllocationToSurvival(0.2),
                   addBinary(true),
                   addAgeSpecific(false),
                   addQuality(true),
                   addInvestmentInRepair(false),
                   addInvestmentAffectingOffspringQuality(true){
                       numOfGametes = maximumAge * numOfOffspringPerFemale;
                   }
    
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
    unsigned int numOfIndividualsToFollow; // number of individuals to follow longitudinal
    unsigned int numOfGametes; // the derived number of gametes a female should have
    float weightMaternalEffect; // to determine how much the maternal quality affects the offspring in ratio to paternal effect.
    float initInvestmentInRepair; // initial investement in repair vs reproduction
    unsigned int numOfStemCells; // number of stem cells a male should create
    float meanMutationBiasInvestmentInRepair; // mean for normal distribution to draw mutation effect on investment in repair genes
    float sdMutationalEffectInvestmentInRepair; // sd for normal distribution to draw mutation effect on investment in repair genes
    float mutationProbInvestmentGenes; // mutation rate of age-specific investment in repair genes
    float weightInvestment; // to weigh the investment in repair genes
    unsigned int scalingParameterForNumOfOffspring; // scaling value for determining number of offspring per individual
    float pointOfHalfScalingParam; // value between 0 and 1 at which half of max offspring is defined
    float baselineSurvival; // baseline survival in the resource distribution when allocation to reproduction is zero
    float scalingStrengthOfAllocationToSurvival; // to determine the strength of the allocative effect on survival
    bool addBinary; // add binary genes to model
    bool addAgeSpecific; // adds age-specific genes to model
    bool addQuality; // adds quality effect to model
    bool addInvestmentInRepair; // adds investment in repair to model
    bool addInvestmentAffectingOffspringQuality; // the investment in repair genes only influence offspring quality
    
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
            checkParam(parID, "meanMutationBias", meanMutationBias, ifs);
            checkParam(parID, "sdMutationalEffectSize", sdMutationalEffectSize, ifs);            
            checkParam(parID, "mutationProbInvestmentGenes", mutationProbInvestmentGenes, ifs);
            checkParam(parID, "sdMutationalEffectInvestmentInRepair", sdMutationalEffectInvestmentInRepair, ifs);
            checkParam(parID, "weightInvestment", weightInvestment, ifs);
            checkParam(parID, "tEnd", tEnd, ifs);
        }
        else break;
    }
}

