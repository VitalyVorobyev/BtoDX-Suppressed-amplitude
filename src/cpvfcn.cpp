#include "cpvfcn.h"

#include <unordered_map>
#include <cmath>
#include <iostream>

#include "mylibs/libTatami/toypdf.h"

#include "cfg.h"
#include "lambda.h"
#include "binnedparams.h"

using std::vector;
using std::abs;
using std::log;
using std::make_pair;
using std::unordered_map;

using std::cout;
using std::endl;

constexpr auto rad_to_deg = 180. / M_PI;

CPVFcn::CPVFcn(fitmode fmode, const vector<BEvt>& evtv, uint16_t sb) :
    mode(fmode), sevtv(evtv), single_bin(sb) {}

double CPVFcn::operator()(const vector<double>& par) const {
    static auto spdf = Cfg::pdf();
    static auto bpars = Cfg::bpars(mode == fitmode::Approx);
    static unordered_map<dtypes, Lambda_b2cud> slambdas;
    if (slambdas.empty()) {  // initialize once
        if (mode != fitmode::Simple) {
            slambdas.emplace(make_pair(dtypes::CPp, Cfg::lamf(dtypes::CPp)));
            slambdas.emplace(make_pair(dtypes::CPn, Cfg::lamf(dtypes::CPn)));
            slambdas.emplace(make_pair(dtypes::KPI, Cfg::lamf(dtypes::KPI)));
            slambdas.emplace(make_pair(dtypes::PIK, Cfg::lamf(dtypes::PIK)));
        } else {
            bpars->setRb(0.);
            slambdas.emplace(make_pair(dtypes::CPp, Cfg::lamf0(dtypes::CPp)));
            slambdas.emplace(make_pair(dtypes::CPn, Cfg::lamf0(dtypes::CPn)));
            slambdas.emplace(make_pair(dtypes::KPI, Cfg::lamf0(dtypes::KPI)));
            slambdas.emplace(make_pair(dtypes::PIK, Cfg::lamf0(dtypes::PIK)));
        }
    }

    for (auto& lamf : slambdas) lamf.second.beta(par[0]);
    bpars->setBeta(par[0]);
    if ((mode == fitmode::Full) || (mode == fitmode::Approx)) {
        for (auto& lamf : slambdas) {
            lamf.second.rb(par[1]);
            lamf.second.delb(par[2]);
        }
        bpars->setRb(par[1]);
        bpars->setDeltaB(par[2]);
    }
    // event loop
    double loglh = 0;
    for (const auto& evt : sevtv) {
        if (single_bin && (abs(evt.Bin()) != single_bin)) continue;
        double ccoef, scoef;
        if (evt.Type() == dtypes::KsPIPI) {
            auto coefs = bpars->coefs(evt.Bin());
            ccoef = coefs.first;
            scoef = coefs.second;
        } else {
            auto &lamf = slambdas.at(evt.Type());
            ccoef = lamf.ccoef();
            scoef = lamf.scoef();
        }
        loglh += log(spdf(evt.t(), evt.Tag(), ccoef, scoef));
    }
    cout << "llh " << loglh << ", beta " << par[0] * rad_to_deg << endl;
    return -2. * loglh;
}

//double CPVFcn::operator()(const vector<double>& par) const {
//    return (par[0]-0.1)*(par[0]-0.1);
//}
