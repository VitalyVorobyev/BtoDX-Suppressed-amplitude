#ifndef LAMBDA_H
#define LAMBDA_H

#include <complex>

class Lambda_b2cud {
 public:
    Lambda_b2cud(double _rd, double _deld, double _rb, double _delb,
                 double _beta, double _gamma);

    std::complex<double> operator()();

    double real() const {return _lambda.real();}
    double imag() const {return _lambda.imag();}

    double ccoef() const {return _ccoef;}
    double scoef() const {return _scoef;}

    void set_rd(double x) {rd = x; update();}
    void set_deld(double x) {deld = x; update();}
    void set_rb(double x) {rb = x; update();}
    void set_delb(double x) {delb = x; update();}
    void set_gamma(double x) {gamma = x; update();}
    void set_beta(double x);

 private:
    void calculate();
    void lambda();
    void update();

    double rd;
    double deld;
    double rb;
    double delb;
    double beta;
    double gamma;

    std::complex<double> _lambda;
    double _ccoef;
    double _scoef;

    static const double deg_to_rag;
};

#endif  // LAMBDA_H
