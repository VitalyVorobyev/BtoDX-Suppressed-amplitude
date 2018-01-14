#include "cpvfitter.h"

#include <fstream>
#include <random> // std::default_random_engine

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"

#include "cfg.h"
#include "cpvfcn.h"
#include "ddfcn.h"

using std::vector;
using std::string;
using std::ifstream;

using std::endl;
using std::cerr;
using std::cout;

using std::pair;

using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnUserParameterState;

constexpr auto rad_to_deg = 180. / M_PI;

CPVFitter::CPVFitter(fitmode mode) : sfmode(mode), single_bin(0) {}

CPVFitter::CPVFitter(string& input_file_name, fitmode mode) : CPVFitter(mode) {
    ReadBEvt(input_file_name);
}

CPVFitter::CPVFitter(vector<string>& input_files, fitmode mode) :
    CPVFitter(mode) {
    for (auto fname : input_files) ReadBEvt(fname);
}

MnUserParameterState CPVFitter::MakeFit() {
    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
    std::default_random_engine rndmeng;
    rndmeng.seed(std::random_device {}());

    MnUserParameters upar;
    upar.Add("beta", udist(rndmeng), 1., -10., 10.);
    if ((sfmode == fitmode::Full) || (sfmode == fitmode::Approx)) {
        upar.Add("rb", 0.02, 0.01, 0., 0.5);
        upar.Add("db", udist(rndmeng), 1., -10., 10.);
    }

    CPVFcn fcn(sfmode, sevtv, single_bin);
    MnMigrad migrad(fcn, upar);
    auto fitmin = migrad();
    auto pstate = fitmin.UserState();

    auto beta = pstate.Value(0) * rad_to_deg;
    while (beta < 0) beta += 180.;
    while (beta > 180) beta -= 180.;
    cout << upar.Name(0) << ": " << upar.Value(0) * rad_to_deg << " -> "
         << beta << " +- " << pstate.Error(0) * rad_to_deg << endl;
    if ((sfmode == fitmode::Full) || (sfmode == fitmode::Approx)) {
        cout << upar.Name(1) << ":   " << upar.Value(1) << " -> "
             << pstate.Value(1) << " +- " << pstate.Error(1) << endl;
        auto db = pstate.Value(2) * rad_to_deg;
        while (db < 0) db += 360.;
        while (db > 360) db -= 360.;
        cout << upar.Name(2) << ":   " << upar.Value(2) * rad_to_deg << " -> "
             << db << " +- " << pstate.Error(2) * rad_to_deg << endl;
    }
    return pstate;
}

auto CPVFitter::ReadDBEvt(const string& fname) {
    vector<DBEvt> evtv;
    ifstream ifile(fname, std::ifstream::in);
    if (!ifile.is_open()) {
        cerr << "Can't open file " << fname << endl;
        return move(evtv);
    }
    while (ifile.good()) {
        static double time;
        static int16_t tag, dbin, bbin;
        ifile >> time >> tag >> dbin >> bbin;
        evtv.emplace_back(time, tag, dbin, bbin);
    }
    cout << evtv.size() << " events found" << endl;
    return move(evtv);
}

MnUserParameterState CPVFitter::DDFit(const string& fname, AbsDDPars& pars) {
    auto evtv = ReadDBEvt(fname);

    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
    std::default_random_engine rndmeng;
    rndmeng.seed(std::random_device {}());

    MnUserParameters upar;
    upar.Add("beta", udist(rndmeng), 1., -10., 10.);

    DDFcn fcn(sfmode, evtv, single_bin, pars);
    MnMigrad migrad(fcn, upar);
    auto fitmin = migrad();
    auto pstate = fitmin.UserState();

    auto beta = pstate.Value(0) * rad_to_deg;
    while (beta < 0) beta += 180.;
    while (beta > 180) beta -= 180.;
    cout << upar.Name(0) << ": " << upar.Value(0) * rad_to_deg << " -> "
         << beta << " +- " << pstate.Error(0) * rad_to_deg << endl;
    return pstate;
}

pair<vector<double>, vector<double>> CPVFitter::Scan(
        double lo, double hi, uint16_t nbin) const {
    CPVFcn fcn(sfmode, sevtv, single_bin);
    pair<vector<double>, vector<double>> r;
    auto dbeta = (hi - lo) / nbin;
    for (auto beta = lo; beta < hi+0.001; beta += dbeta) {
        r.first.emplace_back(beta);
        const vector<double> par{beta / rad_to_deg};
        r.second.emplace_back(fcn(par));
    }
    return r;
}

void CPVFitter::ReadBEvt(string &input_file_name) {
    sevtv.clear();
    ifstream ifile(input_file_name, std::ifstream::in);
    if (!ifile.is_open()) {
        cerr << "Can't open file " << input_file_name << endl;
        return;
    }
    double time;
    int16_t tag;
    int type;
    int16_t bin;
    while (ifile.good()) {
        ifile >> time >> tag >> bin >> type;
        sevtv.push_back(BEvt(time, tag, bin, static_cast<dtypes>(type)));
    }
    ifile.close();
}











