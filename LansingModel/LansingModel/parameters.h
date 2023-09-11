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
    Parameters() : populationSize(1000), // N in manuscript
                   initDamageProportion(0.1),
                   numOfOffspringPerFemale(1), // o in manuscript
                   mutationProb(0.0024), // mu_b in manuscript
                   extrinsicMortRisk(0.0),
                   outputTime(10),
                   tEnd(10000),
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
                   includeRecombination(false),
                   addBinary(false),
                   addAgeSpecific(false),
                   addQuality(true),
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
    bool includeRecombination; // bool to either add recombination or not
    bool addBinary; // add binary genes to model
    bool addAgeSpecific; // adds age-specific genes to model
    bool addQuality; // adds quality effect to model
    bool addInvestmentInRepair; // adds investment in repair to model
    std::string temp_params_to_record; // to temporarily store varying parameter(s)
    std::vector < std::string > param_names_to_record; // to store the names of the varying parameter(s)
    std::vector < float > params_to_record; // to store the varying parameter(s)
    
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
    void checkParam(const std::string parID,
                    const std::string focal_parametername,
                    std::string& parameter,
                    std::ifstream& ifs);
    void setAdditionalParams();
    std::vector<std::string> split(std::string s);
    std::vector<float> create_params_to_record(const std::vector< std::string >& param_names);
    float get_val(std::string s);
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

void Parameters::checkParam(const std::string parID,
                            const std::string focal_parametername,
                            std::string& parameter,
                            std::ifstream& ifs) {
    // set parameter from file to parameter in object bool parameter
    if (parID == focal_parametername) {
        ifs >> parameter;
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
            checkParam(parID, "params_to_record", temp_params_to_record, ifs);
            checkParam(parID, "includeRecombination", includeRecombination, ifs);
            checkParam(parID, "scalingStrengthOfAllocationToReproduce", scalingStrengthOfAllocationToReproduce, ifs);
            checkParam(parID, "steepnessAllocationToReproduce", steepnessAllocationToReproduce, ifs);
            
        }
        else break;
    }
    
    param_names_to_record = split(temp_params_to_record);
    params_to_record = create_params_to_record(param_names_to_record);
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

std::vector< std::string > Parameters::split(std::string s) {
    // code from: https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
    std::vector< std::string > output;
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      output.push_back(token);
      s.erase(0, pos + delimiter.length());
    }
    output.push_back(s);
    return output;
}

std::vector< float > Parameters::create_params_to_record(const std::vector< std::string >& param_names) {
  std::vector< float > output;
  for (auto i : param_names) {
    output.push_back(get_val(i));
  }
  return output;
}


float Parameters::get_val(std::string s) {
    if (s == "addBinary")                               return addBinary;
    if (s == "addQuality")                              return addQuality;
    if (s == "addAgeSpecific")                          return addAgeSpecific;
    if (s == "addInvestmentInRepair")                   return addInvestmentInRepair;
    if (s == "populationSize")                          return populationSize;
    if (s == "tEnd")                                    return tEnd;
    if (s == "mutationProb")                            return mutationProb;
    if (s == "mutationProbStemcell")                    return mutationProbStemcell;
    if (s == "mutationProbAgeSpecificGenes")            return mutationProbAgeSpecificGenes;
    if (s == "meanMutationBias")                        return meanMutationBias;
    if (s == "sdMutationalEffectSize")                  return sdMutationalEffectSize;
    if (s == "mutationProbInvestmentGenes")             return mutationProbInvestmentGenes;
    if (s == "includeRecombination")                    return includeRecombination;
    if (s == "weightInvestment")                        return weightInvestment;
    if (s == "scalingStrengthOfAllocationToReproduce")  return scalingStrengthOfAllocationToReproduce;
    if (s == "steepnessAllocationToReproduce")          return steepnessAllocationToReproduce;
    

    throw std::runtime_error("can not find parameter");
    return -1.f; // FAIL
}
