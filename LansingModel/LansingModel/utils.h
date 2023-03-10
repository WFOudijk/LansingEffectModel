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

#endif // UTILS_H
