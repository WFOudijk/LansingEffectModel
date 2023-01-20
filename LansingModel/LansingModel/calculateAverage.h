//
//  calculateAverage.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 17/01/2023.
//

#pragma once
//
//using indVec = std::vector<Individual>;
//using arrayOfGenes = std::array<double, numAgeClasses>;
//
//arrayOfGenes calcAverageAcrossAgeClasses(const indVec& individuals) {
//    /**Function to calculate the average survival probability across the age classes
//     over a vector of individuals. */
//    arrayOfGenes popAverage;
//    popAverage.fill(0.0);
//
//    for (const auto& i : individuals) { // loop through the individuals
//        for(auto j = 0u; j < i.averageSurvivalProb.size(); ++j){ // loop through every age group
//            // sum the survival prob average of every individual for the jth age group
//            popAverage[j] += i.averageSurvivalProb[j];
//        }
//    }
//    for(auto& i : popAverage) { // loop through the averages
//        i *= 1.0 / (individuals.size()); // compare the sum per age group with the size of the individuals
//    }
//    return popAverage; // return array with average survival probability for every age over the individuals
//}
//
//std::vector<double> calcLifeExpectancyPerIndividual(const indVec& individuals){
//    /**Function to calculate the Life Expectancy per individual . **/
//    std::vector<double> lifeExpectancy; // for actual life expectancy per individual
//    for (auto individual = 0u; individual < individuals.size(); ++individual){ // loop through every individual
//        std::vector<double> lifeExpectancyPerIndividual; // per individual, a likelihood for every age class
//        for (auto i = 0u; i < individuals[individual].genesMaternal.size(); ++i){ // loop through every age
//            double lifeExpectancyInd = // get survival prob of current age and individual
//                    individuals[individual].averageSurvivalProb[i];
//            for (int j = 1; j <= i; ++j){ // to make sure every previous survival prob is taken into account
//                // current survivalprob multiplied with every previous survivalprob (per age) of this individual
//                lifeExpectancyInd *=
//                    individuals[individual].averageSurvivalProb[i - j];
//            }
//            lifeExpectancyPerIndividual.push_back(lifeExpectancyInd); // add the product to a vector
//        }
//        // the vector is the numOfAge long, with one value per age class. This value represents the likelihood
//        // of the individual becoming that specific age. By summing this you get the life expectancy.
//        double sum = std::accumulate(lifeExpectancyPerIndividual.begin(),
//                                     lifeExpectancyPerIndividual.end(), 0.0);
//        lifeExpectancy.push_back(sum); // this sum is added to a vector
//    }
//    return lifeExpectancy; // this vector consists of the life expectancy per individual
//}
