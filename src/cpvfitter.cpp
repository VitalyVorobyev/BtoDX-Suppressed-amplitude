#include "cpvfitter.h"

#include <fstream>
#include <random> // std::default_random_engine

#include "Minuit2/MnUserParameters.h"
#include "Minuit2/MnMigrad.h"
#include "Minuit2/FunctionMinimum.h"
#include "Minuit2/MnStrategy.h"
//#include "Minuit2/MnMachinePrecision.h"

#include "Math/MinimizerOptions.h"

#include "cfg.h"
#include "cpvfcn.h"

using std::vector;
using std::string;
using std::ifstream;

using std::endl;
using std::cerr;
using std::cout;

using ROOT::Minuit2::MnMigrad;
using ROOT::Minuit2::MnUserParameters;
using ROOT::Minuit2::MnUserParameterState;
using ROOT::Minuit2::MnStrategy;

constexpr auto rad_to_deg = 180. / M_PI;

CPVFitter::CPVFitter(fitmode mode) : sfmode(mode) {}

CPVFitter::CPVFitter(string& input_file_name, fitmode mode) : CPVFitter(mode) {
    ReadData(input_file_name);
}

CPVFitter::CPVFitter(vector<string>& input_files, fitmode mode) :
    CPVFitter(mode) {
    for (auto fname : input_files) ReadData(fname);
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

    MnStrategy mnstr(2);
//    mnstr.setHighStrategy();

    CPVFcn fcn(sfmode, sevtv, single_bin);
    MnMigrad migrad(fcn, upar);
    auto fitmin = migrad();
    auto pstate = fitmin.UserState();

//    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1.e-6);

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

void CPVFitter::ReadData(string &input_file_name) {
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
