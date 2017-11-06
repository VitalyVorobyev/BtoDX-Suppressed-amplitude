#ifndef RNDMD_H
#define RNDMD_H

#include <random>  // std::default_random_engine

class RndmD {
    static std::default_random_engine rndmeng;
    const double m_lo;
    const double m_hi;
    std::uniform_real_distribution<double> udist;

 public:
    RndmD(double lo, double hi) :
        m_lo(lo), m_hi(hi), udist(lo, hi) {}
    double operator()() {return udist(rndmeng);}
    static void setSeed(int seed=0) {
        if (seed) rndmeng.seed(seed);
        else rndmeng.seed(std::random_device {}());
    }
};

#endif // RNDMD_H
