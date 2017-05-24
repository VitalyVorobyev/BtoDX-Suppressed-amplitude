#include "cpvfitter.h"

#include <fstream>
#include <cmath>
#include <random> // std::default_random_engine

#include "TVirtualFitter.h"

typedef TVirtualFitter VFit;

using std::endl;
using std::cerr;
using std::cout;
using std::log;

const double rad_to_deg = 180. / M_PI;

std::unordered_map<dtypes, Lambda_b2cud>* CPVFitter::slambdas = nullptr;
libTatami::ToyPdf* CPVFitter::spdf = nullptr;
std::vector<Evt>* CPVFitter::sevtv = nullptr;
fitmode CPVFitter::sfmode = fitmode::Simple;

CPVFitter::CPVFitter(fitmode mode) :
    fmode(mode), pdf(Cfg::pdf()) {
    sfmode = fmode;
    spdf = &pdf;
    if (mode != fitmode::Simple) {
        lambdas.insert({dtypes::CPp, Cfg::lamf(dtypes::CPp)});
        lambdas.insert({dtypes::CPn, Cfg::lamf(dtypes::CPn)});
        lambdas.insert({dtypes::KPI, Cfg::lamf(dtypes::KPI)});
        lambdas.insert({dtypes::PIK, Cfg::lamf(dtypes::PIK)});
    } else {
        lambdas.insert({dtypes::CPp, Cfg::lamf0(dtypes::CPp)});
        lambdas.insert({dtypes::CPn, Cfg::lamf0(dtypes::CPn)});
        lambdas.insert({dtypes::KPI, Cfg::lamf0(dtypes::KPI)});
        lambdas.insert({dtypes::PIK, Cfg::lamf0(dtypes::PIK)});
    }
    slambdas = &lambdas;
}

CPVFitter::CPVFitter(std::string& input_file_name, fitmode mode) :
    CPVFitter(mode) {
    ReadData(input_file_name);
    if (sevtv->size() > 0) MakeFit();
}

CPVFitter::CPVFitter(std::vector<std::string>& input_files, fitmode mode) :
    CPVFitter(mode) {
    for (auto fname : input_files) ReadData(fname);
    if (sevtv->size() > 0) MakeFit();
}

void CPVFitter::TheFCN(int& npar, double* grad,
                       double &fval, double *p, int iflag) {
    for (auto& lamf : *slambdas) {
        lamf.second.set_beta(*p);
        if (sfmode == fitmode::Full) {
            lamf.second.set_rb(*(p+1));
            lamf.second.set_delb(*(p+2));
        }
    }
    double loglh = 0;
    for (auto& evt : *sevtv) {
        auto &lamf = slambdas->find(evt.type)->second;
        loglh += -2.*log((*spdf)(evt.time, evt.tag,
                                 lamf.ccoef(), lamf.scoef()));
    }
    cout << "llh " << loglh << ", beta " << (*p) * rad_to_deg << endl;
    fval = loglh;
}

void CPVFitter::MakeFit() {
    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
    std::default_random_engine rndmeng;
    rndmeng.seed(std::random_device {}());

    VFit::SetDefaultFitter("Minuit");
    const int npar = fmode == fitmode::Full ? 3 : 1;
    VFit* fitter = VFit::Fitter(0, npar);
    fitter->SetParameter(0, "beta", udist(rndmeng), 1., 0., 0.);
    if (fmode == fitmode::Full) {
        fitter->SetParameter(1, "rb", 0.02, 0.01, 0., .2);
        fitter->SetParameter(2, "db", udist(rndmeng), 1., -2.*M_PI, 2.*M_PI);
    }
    fitter->SetFCN(TheFCN);

    double arglist[100];
    // Set pring level
    arglist[0] = 0;
    fitter->ExecuteCommand("SET PRINT", arglist, 2);
    // Minimize
    arglist[0] = 0;  // number of calls
    arglist[1] = 0.01;  // tolerance
    fitter->ExecuteCommand("MIGRAD", arglist, 2);
    // Get results
    double beta_fit, beta_error;
    beta_fit = fitter->GetParameter(0) * rad_to_deg;
    if (beta_fit < 0) beta_fit += 180.;
    while (beta_fit > 90) beta_fit -= 90.;
    beta_error = fitter->GetParError(0) * rad_to_deg;
    cout << "beta = " << beta_fit << " +- " << beta_error << endl;

    double rb_fit = fitter->GetParameter(1);
    double rb_error = fitter->GetParError(1);
    cout << "  rb = " << rb_fit << " +- " << rb_error << endl;

    double db_fit = fitter->GetParameter(2) * rad_to_deg;
    double db_error = fitter->GetParError(2) * rad_to_deg;
    cout << "delb = " << db_fit << " +- " << db_error << endl;

//    double chi2, edm, errdef;
//    int nvpar, nparx;
//    fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
//    cout << " loglh: " << chi2 << endl
//         << "   edm: " << edm << endl
//         << "errded: " << errdef << endl
//         << " nvpar: " << nvpar << endl
//         << " nparx: " << nparx << endl;
}

void CPVFitter::ReadData(std::string &input_file_name) {
    std::ifstream ifile;
    ifile.open(input_file_name, std::ifstream::in);
    if (!ifile.is_open()) {
        cerr << "Can't open file " << input_file_name << endl;
        return;
    }
    double time;
    int tag;
    int type;
    while (ifile.good()) {
        ifile >> time >> tag >> type;
        evtv.push_back(Evt(time, tag, static_cast<dtypes>(type)));
    }
    ifile.close();
    sevtv = &evtv;
}
