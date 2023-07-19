//
//  randomizer.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//

#pragma once

struct Randomizer {
    std::mt19937_64 rng;
    Randomizer() {
        unif = std::uniform_real_distribution<double>(0, 1);
    }
    double uniform() {
        return unif(rng);
    }
    void setSeed(size_t s){
        rng.seed(s);
    }
    
    int drawRandomNumber(size_t sizeVector) {
        // to pick a random index in a vector
        if (sizeVector < 1) {
            return 0;
        }
        std::uniform_int_distribution<int> d(0, static_cast<int>(sizeVector) - 1);
        return d(rng);
    }
        
    bool bernoulli(double p = 0.5) {
        return std::bernoulli_distribution(p)(rng);        
    }
    
    double drawMutationEffect() {
            return mutationEffect(rng);
        }
    
    void setMutationEffect(double m, double sd) {
        // create normal distribution based on user defined mean and standard deviation
        mutationEffect = std::normal_distribution<double>(m, sd);
    }
				
    int drawNumOfMuts() {
        return distMutEvents(rng);
    }
				
    void setDistMutEvents(double setDist) {
        distMutEvents = std::poisson_distribution<int>(setDist);
    }
    
    double drawMutationEffectInvestment(){
        return mutationEffectInvestment(rng);
    }
    
    void setMutationEffectInvestment(double m, double sd){
        mutationEffectInvestment = std::normal_distribution<double>(m, sd);
    }
				
    template <typename T>
    T rn(const T n){
        static std::uniform_int_distribution<T> d{};
        using parm_t = typename decltype(d)::param_type;
        return d(rng, parm_t{0,n-1});
    }

    unsigned rui32(){
        static std::uniform_int_distribution<unsigned> d{};
        using parm_t = typename decltype(d)::param_type;
        return d(rng, parm_t{0,UINT32_MAX});
    }

    unsigned long rui64(){
        static std::uniform_int_distribution<unsigned long> d{};
        using parm_t = typename decltype(d)::param_type;
        return d(rng, parm_t{0,UINT64_MAX});
    }

    unsigned rpois(const double lambda){
        static std::poisson_distribution<unsigned> d{};
        using parm_t = typename decltype(d)::param_type;
        return d(rng, parm_t{lambda});
    }
    
     
    std::uniform_real_distribution<double> unif;
    std::normal_distribution<double> mutationEffect;
    std::poisson_distribution<int> distMutEvents;
    std::normal_distribution<double> mutationEffectInvestment;
};
