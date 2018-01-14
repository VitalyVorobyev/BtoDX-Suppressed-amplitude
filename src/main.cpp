#include <utility>
#include <iostream>
#include <string>
#include <map>

#include "wfdriver.h"
#include "cfg.h"

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::map;

const map<string, fitmode> mode_map = {
    {"corr", fitmode::Corrected},
    {"full", fitmode::Full},
    {"simp", fitmode::Simple},
    {"appr", fitmode::Approx}
};

const map<string, double> pars = {
    {"x", 0.01},
    {"y", 0.01},
    {"rb", 0.2},
    {"beta", 23.},
    {"dtlim", 70.}
};

int main(int argc, char** argv) {
    const auto x = pars.at("x");
    const auto y = pars.at("y");

    Cfg::set_beta(pars.at("beta"));
    Cfg::set_dtlim(pars.at("dtlim"));
    Cfg::set_charm_mix(x, y);
    Cfg::print_config();

    const auto bdecay = bmode::DKs;
    const double delb = 0.;
    constexpr uint64_t nevt = 100000;
    Cfg::set_delb(delb);
    cout << "New delb " << delb << endl;
    WFDriver::gen_dd_with_wf(bdecay, pars.at("rb"), nevt, 0, 100);
    WFDriver::gen_cp_with_wf(bdecay, pars.at("rb"), 0.5 * nevt, 0.5 * nevt, 0, 100);
    WFDriver::fit_with_wf(bdecay, fitmode::DKs, pars.at("rb"), 0, 100);

//    for (double delb = 0.; delb < 360.; delb += 10) {
//        Cfg::set_delb(delb);
//        cout << "New delb " << delb << endl;
//        WFDriver::gen_dd_with_wf(bdecay, pars.at("rb"), nevt, 0, 100);
//        WFDriver::gen_cp_with_wf(bdecay, pars.at("rb"), 0.5 * nevt, 0.5 * nevt, 0, 100);
//        WFDriver::fit_with_wf(bdecay, fitmode::DKs, pars.at("rb"), 0, 100);
//    }

//    constexpr uint16_t sb = 0;
//    for (double delb = 0.; delb < 360.; delb += 10) {
//        Cfg::set_delb(delb);
//        cout << "New delb " << delb << endl;
//        WFDriver::gen_dd_with_wf(bdecay, pars.at("rb"), nevt);
//        WFDriver::gen_cp_with_wf(bdecay, pars.at("rb"), 100, nevt);
//        WFDriver::fit_with_wf(bdecay, fitmode::Corrected, pars.at("rb"));
//        WFDriver::fit_with_wf(bdecay, fitmode::Full, pars.at("rb"));
//    }

//   WFDriver::gen_dd_with_charm_mix(bdecay, pars.at("x"), pars.at("y"), nevt);
//   WFDriver::gen_cp_with_charm_mix(bdecay, pars.at("x"), pars.at("y"), nevt, 100);
//   WFDriver::fit_with_charm_mix(bdecay, fitmode::Corrected, pars.at("x"), pars.at("y"));
//   WFDriver::fit_with_charm_mix(bdecay, fitmode::Full, pars.at("x"), pars.at("y"));

//    constexpr uint16_t idxlo = 0;
//    constexpr uint16_t idxhi = 100;

//    double rb = 0.02;
//    WFDriver::gen_cp_with_charm_mix(10, 10, x, y);
//    WFDriver::gen_dd_with_charm_mix(nevt, x, y);
//    WFDriver::fit_with_charm_mix(fitmode::Corrected, x, y, sb);
//    WFDriver::fit_with_charm_mix(fitmode::Full, x, y, sb);

//    WFDriver::gen_dd_with_wf(nevt, idxlo, idxhi, rb);
//    WFDriver::gen_cp_with_wf(0, nevt, idxlo, idxhi, rb);
//    WFDriver::fit_with_wf(fitmode::Corrected, idxlo, idxhi, rb, sb);
//    WFDriver::fit_with_wf(fitmode::Full, idxlo, idxhi, rb, sb);

    return 0;
}
