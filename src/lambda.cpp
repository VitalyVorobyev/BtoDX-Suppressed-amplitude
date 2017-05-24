#include "lambda.h"

#include <cmath>
#include <iostream>

typedef std::complex<double> compld;
typedef Lambda_b2cud Lambda;

using std::exp;
using std::pow;
using std::cout;
using std::endl;

const double Lambda::deg_to_rag = M_PI / 180.;

Lambda::Lambda_b2cud(double _rd, double _deld, double _rb, double _delb,
               double _beta, double _gamma) :
    rd(_rd), deld(deg_to_rag*_deld),
    rb(_rb), delb(deg_to_rag*_delb),
    beta(deg_to_rag*_beta), gamma(deg_to_rag*_gamma) {
    update();
}

void Lambda::lambda() {
    _lambda = rd * exp(compld(0., -deld)) +
              rb * (             exp(compld(0., delb - gamma)) -
                    pow(rd, 2) * exp(compld(0., delb + gamma - 2 * deld))) -
              rd * pow(rb, 2) * exp(compld(0., 2. * (delb - deld)));
    _lambda *= exp(compld(0., -2.*beta));
}

void Lambda::calculate() {
    const double lamfsq = pow(_lambda.imag(), 2) + pow(_lambda.real(), 2);
    _ccoef = (1. - lamfsq) / (1. + lamfsq );
    _scoef = 2. * _lambda.imag() / (1. + lamfsq );
}

void Lambda::update() {
    lambda(); calculate();
}

void Lambda::set_beta(double x) {
    _lambda *= exp(compld(0., 2.*beta));
    beta = x;
    _lambda *= exp(compld(0., -2.*beta));
    calculate();
}
