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
    if (x<0.0) x=0.0; else if (x>1.0) x=1.0;
}

template<typename T>
inline T sqr(T x){
    return x*x;
}

int stochasticRound(const double x, Randomizer& rng) {
    double intpart;
    double fractpart = modf(x, &intpart);
    if (rng.bernoulli(fractpart)) return static_cast<int>(x + 1.0);
    return static_cast<int>(intpart);
}

float logistic(double k, double x, double i){
    return (1.0 / (1.0 + exp(-k * x - i)));
}

#endif // UTILS_H
