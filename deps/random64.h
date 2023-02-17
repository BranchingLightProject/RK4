/**
 * Random number generator
 * Author : Jose Daniel Munoz et al.
 */
#ifndef __RANDOM64_H
#define __RANDOM64_H
#include <cmath>

class CRandom {
    unsigned long long u, v, w;

   public:
    CRandom(unsigned long long seed);
    unsigned long long int64(void);
    unsigned int int32() { return (unsigned int)int64(); };
    double exponencial(float tau);
    double gauss(float mu, float sigma);

    double r(void) { return 5.42101086242752217E-20 * int64(); };
};

#endif  // __RANDOM64_H