#include "abslambda.h"

//#include <cmath>

using compld = std::complex<double>;

using std::pow;
using std::exp;

//const double AbsLambda::deg_to_rag = M_PI / 180.;

void AbsLambda::lambda() {
    m_lambda = afbar() / af() * exp(compld(0., -2.*m_beta));
}

void AbsLambda::calculate() {
    const double lamfsq = pow(m_lambda.imag(), 2) + pow(m_lambda.real(), 2);
    m_ccoef = (1. - lamfsq) / (1. + lamfsq);
    m_scoef = 2. * m_lambda.imag() / (1. + lamfsq);
}

void AbsLambda::beta(double x) {
    m_lambda *= exp(compld(0., 2.*m_beta));
    m_beta = x;
    m_lambda *= exp(compld(0., -2.*m_beta));
    calculate();
}
