#include "driver.h"

#include <iostream>
#include <string>
#include <utility>  // pair
#include <cmath>

void Driver::delb_scan(dtypes type, fitmode mode) {
    std::cout << "delb_scan" << std::endl;
    const int N = pow(10, 6);
    std::string fname = Cfg::dfile(type);
    CPVGen gen;
    for (double db = 0.; db < 360.; db += 15.) {
        std::cout << "db " << db << std::endl;
        Cfg::set_delb(db);
        gen.GenAndWrite(type, std::make_pair(N / 2, N / 2));
        CPVFitter fitter(fname, mode);
    }
}

void Driver::single_fit(dtypes type, fitmode mode) {
    auto lambda = Cfg::lamf(type);
    auto pdf = Cfg::pdf();
    pdf.SetS(lambda.scoef());
    pdf.SetC(lambda.ccoef());
    const int N = pow(10, 6);
    CPVGen gen;
    auto fname = Cfg::dfile(type);
    gen.GenAndWrite(pdf, std::make_pair(N / 2, N / 2), fname, type);
    CPVFitter fitter(fname, mode);
}

void Driver::cpp_and_cpn(fitmode mode) {
    CPVGen gen;
    const int N = pow(10, 6);
    gen.GenAndWrite(dtypes::CPp, std::make_pair(N / 2, N / 2));
    gen.GenAndWrite(dtypes::CPn, std::make_pair(N / 2, N / 2));
    std::vector<std::string> files = {
        Cfg::dfile(dtypes::CPp),
        Cfg::dfile(dtypes::CPn)
    };
    CPVFitter fitter(files, mode);
}

void Driver::all_in_one(fitmode mode) {
    Cfg::set_ncp(2*std::pow(10, 6));
    Cfg::set_nfl(2*std::pow(10, 6));
    CPVGen gen;
    gen.CompleteData();
    std::vector<std::string> files = {
        Cfg::dfile(dtypes::CPp),
        Cfg::dfile(dtypes::CPn),
        Cfg::dfile(dtypes::KPI),
        Cfg::dfile(dtypes::PIK)
    };
    CPVFitter fitter(files, mode);
}
