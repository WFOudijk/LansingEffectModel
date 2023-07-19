//
//  utils.h
//  LansingModel
//
//  Created by Willemijn Oudijk on 24/02/2023.
//

#ifndef UTILS_H
#define UTILS_H


template<typename T=double>
inline void clip01(T& x){
    /**Function to clip the input if it decreases 0.0 or exceeds 1.0 to the respective limit**/
    if (x<0.0) x=0.0; else if (x>1.0) x=1.0;
}

template<typename T>
inline T sqr(T x){
    /**Function to square the input **/
    
    return x*x;
}

float logistic(double k, double x, double i){
    /**Function to use when a logistic function is needed. **/
    
    return (1.0 / (1.0 + exp(-k * x - i)));
}

#endif // UTILS_H
