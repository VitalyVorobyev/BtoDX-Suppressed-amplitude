#ifndef LAMBDA_H
#define LAMBDA_H

#include "abslambda.h"

#include <complex>

class Lambda_b2cud : public AbsLambda {
    std::complex<double> af() const override;
    std::complex<double> afbar() const override;
    void update() {lambda(); calculate();}

    double m_rd;
    double m_deld;
    double m_rb;
    double m_delb;
    double m_gamma;

 public:
    Lambda_b2cud(double _rd, double _deld, double _rb, double _delb,
                 double _beta, double _gamma);
    void rd(double x) {m_rd = x; update();}
    void deld(double x) {m_deld = x; update();}
    void rb(double x) {m_rb = x; update();}
    void delb(double x) {m_delb = x; update();}
    void gamma(double x) {m_gamma = x; update();}
};

#endif  // LAMBDA_H
