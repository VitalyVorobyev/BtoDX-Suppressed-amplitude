#include <utility>
#include <iostream>

#include "driver.h"
#include "cfg.h"
#include "cpvgen.h"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
    Cfg::print_config();
//    Driver driver;
//    driver.delb_scan(dtypes::CPn);
//    driver.single_fit(dtypes::CPp);
//    driver.single_fit(dtypes::CPn);
//    driver.single_fit(dtypes::KPI, fitmode::Corrected);
//    driver.single_fit(dtypes::PIK);
//    driver.cpp_and_cpn(fitmode::Simple);
//    driver.cpp_and_cpn(fitmode::Corrected);
//    driver.cpp_and_cpn(fitmode::Full);
//    for (unsigned i = 0; i < 500; i++)
//        driver.all_in_one(fitmode::Corrected);
//    driver.all_in_one(fitmode::Simple);

//    CPVGen gen;
//    auto nevt = std::make_pair(1000000, 1000000);
//    gen.GenAndWrite(dtypes::KsPIPI, nevt);

//    auto fname = Cfg::dfile(dtypes::KsPIPI);
//    CPVFitter fit(fname, fitmode::Full);
//    CPVFitter fit(fname, fitmode::Corrected);
//    CPVFitter fit1(fname, fitmode::Simple);
//    CPVFitter fit2(fname, fitmode::Approx);
    auto nevt = 100000;
    CPVGen gen;
    Driver driver;
    auto nevts = std::make_pair(nevt, nevt);

    // Null scan
    Cfg::set_rb(0.);
    cout << "Null scan: rB = 0" << endl;
    gen.GenAndWrite(dtypes::KsPIPI, nevts);

    cout << "### FULL DP SIMPLE ###" << endl;
    driver.kspp_fit(fitmode::Simple, nevt);

//    for(auto idx = 0; idx < 10; idx++) {
//        cout << "Null scan: rB = 0" << endl;
//        gen.GenAndWrite(dtypes::KsPIPI, nevts);

//        cout << "### BIN SCAN SIMPLE ###" << endl;
//        driver.kspp_bin_scan(fitmode::Simple, nevt);
//    }

    // Simple scan
//    for (auto delb = 0; delb < 360; delb += 15) {
//        cout << "### DELTAB " << delb << endl;
//        Cfg::set_delb(delb);
//        gen.GenAndWrite(dtypes::KsPIPI, nevts);

//        cout << "### BIN SCAN SIMPLE ###" << endl;
//        driver.kspp_bin_scan(fitmode::Simple, nevt);

//        cout << "### FULL DP SIMPLE ###" << endl;
//        driver.kspp_fit(fitmode::Simple, nevt);
//    }

    // Approx scan
//    for (auto delb = 0; delb < 360; delb += 15) {
//        cout << "### DELTAB " << delb << endl;
//        Cfg::set_delb(delb);
//        gen.GenAndWrite(dtypes::KsPIPI, nevts);

//        cout << "### FULL DP APPROX ###" << endl;
//        driver.kspp_fit(fitmode::Approx, nevt);
//    }

    // Full scan
//    for (auto delb = 0; delb < 360; delb += 15) {
//        cout << "### DELTAB " << delb << " ###" << endl;
//        Cfg::set_delb(delb);
//        gen.GenAndWrite(dtypes::KsPIPI, nevts);

//        cout << "### FULL DP APPROX ###" << endl;
//        driver.kspp_fit(fitmode::Full, nevt);
//    }
    return 0;
}
