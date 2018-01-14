#include "lambda.h"

#include <cmath>

using compld = std::complex<double>;
using Lambda = Lambda_b2cud;

using std::exp;

Lambda::Lambda_b2cud(double _rd, double _deld, double _rb, double _delb,
               double _beta, double _gamma) : AbsLambda(_beta),
    m_rd(_rd), m_deld(deg_to_rag * _deld),
    m_rb(_rb), m_delb(deg_to_rag * _delb),
    m_gamma(deg_to_rag * _gamma) {
    update();
}

compld Lambda::af() const {
    return 1. + m_rd * m_rb * exp(compld(0., -m_deld + m_delb + m_gamma));
}

compld Lambda::afbar() const {
    return m_rd * exp(compld(0., -m_deld)) +
           m_rb * exp(compld(0.,  m_delb - m_gamma));
}

//void Lambda::lambda() {
//    _lambda = rd * exp(compld(0., -deld)) +
//              rb * (             exp(compld(0., delb - gamma)) -
//                    pow(rd, 2) * exp(compld(0., delb + gamma - 2 * deld))) -
//              rd * pow(rb, 2) * exp(compld(0., 2. * (delb - deld)));
//    _lambda *= exp(compld(0., -2.*beta));
//}
