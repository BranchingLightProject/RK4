#include"random64.h"

CRandom::CRandom(unsigned long long seed){
    v = 4101842887655102017LL; w=1;
    u = seed ^ v; int64();
    v = u; int64();
    w = v; int64();
}
unsigned long long CRandom::int64(void){
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; v ^= v << 31; v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long x = u ^ (u << 21); x ^= x >> 35; x ^= x << 4;
    
    return (x + v) ^ w;
}
inline double CRandom::exponencial(float tau){
    return -tau*std::log(r());
}
inline double CRandom:: gauss(float mu, float sigma){
    return sigma*std::sqrt(-2*std::log(r()))*std::cos(2*M_PI*r()) + mu;
}


