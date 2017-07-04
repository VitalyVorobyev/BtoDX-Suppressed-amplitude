#ifndef ABSLAMBDA_H
#define ABSLAMBDA_H

#include <complex>
#include <cmath>

class AbsLambda {
    virtual std::complex<double> af() const = 0;
    virtual std::complex<double> afbar() const = 0;

    std::complex<double> m_lambda;
    double m_ccoef;
    double m_scoef;

    double m_beta;

 protected:
    void calculate();
    void lambda();
    static constexpr auto deg_to_rag = M_PI / 180.;

public:
    AbsLambda(double _beta) : m_beta(deg_to_rag*_beta) {}

    double beta() const {return m_beta;}
    void beta(double x);

    std::complex<double> operator()();

    auto real() const {return m_lambda.real();}
    auto imag() const {return m_lambda.imag();}

    auto ccoef() const {return m_ccoef;}
    auto scoef() const {return m_scoef;}
};

#endif  // ABSLAMBDA_H
