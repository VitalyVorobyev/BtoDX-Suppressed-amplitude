#include "minuitfitter.h"

//#include <fstream>
//#include <random> // std::default_random_engine
//#include <unordered_map>
//#include <cmath>
//#include <map>
//#include <iostream>

//#include "cfg.h"
//#include "lambda.h"
//#include "binnedparams.h"

//#include "TVirtualFitter.h"

//using VFit = TVirtualFitter;
//using DDMap = std::unordered_map<int16_t,
//              std::unordered_map<int16_t, std::pair<double, double>>>;

//using std::vector;
//using std::string;
//using std::make_pair;
//using std::unordered_map;

//using std::ifstream;

//using std::endl;
//using std::cerr;
//using std::cout;

//using std::log;
//using std::abs;

//using std::unique_ptr;

//constexpr auto rad_to_deg = 180. / M_PI;

//vector<BEvt> MinuitFitter::sevtv;
//vector<DBEvt> MinuitFitter::sdevtv;
//uint16_t MinuitFitter::single_bin;
//fitmode MinuitFitter::sfmode;
//std::unique_ptr<AbsDDPars> MinuitFitter::m_pars;
//AbsDDPars MinuitFitter::m_pars;

//MinuitFitter::MinuitFitter(fitmode mode) {
//    sfmode = mode;
//}

//MinuitFitter::MinuitFitter(const string& input_file_name, fitmode mode) :
//    MinuitFitter(mode) {
//    ReadData(input_file_name);
//}

//MinuitFitter::MinuitFitter(const vector<string>& input_files, fitmode mode) :
//    MinuitFitter(mode) {
//    for (const auto& fname : input_files) ReadData(fname);
//}

//void MinuitFitter::TheFCN(int& npar, double* grad,
//                          double &fval, double *p, int iflag) {
//    static auto spdf = Cfg::pdf();
//    static auto bpars = Cfg::bpars(sfmode == fitmode::Approx);
//    static unordered_map<dtypes, Lambda_b2cud> slambdas;
//    if (slambdas.empty()) {  // initialize once
//        if (sfmode != fitmode::Simple) {
//            slambdas.emplace(make_pair(dtypes::CPp, Cfg::lamf(dtypes::CPp)));
//            slambdas.emplace(make_pair(dtypes::CPn, Cfg::lamf(dtypes::CPn)));
//            slambdas.emplace(make_pair(dtypes::KPI, Cfg::lamf(dtypes::KPI)));
//            slambdas.emplace(make_pair(dtypes::PIK, Cfg::lamf(dtypes::PIK)));
//        } else {
//            bpars->setRb(0.);
//            slambdas.emplace(make_pair(dtypes::CPp, Cfg::lamf0(dtypes::CPp)));
//            slambdas.emplace(make_pair(dtypes::CPn, Cfg::lamf0(dtypes::CPn)));
//            slambdas.emplace(make_pair(dtypes::KPI, Cfg::lamf0(dtypes::KPI)));
//            slambdas.emplace(make_pair(dtypes::PIK, Cfg::lamf0(dtypes::PIK)));
//        }
//    }

//    for (auto& lamf : slambdas) lamf.second.beta(*p);
//    bpars->setBeta(*p);
//    if ((sfmode == fitmode::Full) || (sfmode == fitmode::Approx)) {
//        for (auto& lamf : slambdas) {
//            lamf.second.rb(*(p+1));
//            lamf.second.delb(*(p+2));
//        }
//        bpars->setRb(*(p+1));
//        bpars->setDeltaB(*(p+2));
//    }
//    // event loop
//    double loglh = 0;
//    for (const auto& evt : sevtv) {
//        if (single_bin && (abs(evt.Bin()) != single_bin)) continue;
//        double ccoef, scoef;
//        if (evt.Type() == dtypes::KsPIPI) {
//            auto coefs = bpars->coefs(evt.Bin());
//            ccoef = coefs.first;
//            scoef = coefs.second;
//        } else {
//            auto &lamf = slambdas.find(evt.Type())->second;
//            ccoef = lamf.ccoef();
//            scoef = lamf.scoef();
//        }
//        loglh += -2. * log(spdf(evt.t(), evt.Tag(), ccoef, scoef));
//    }
////    cout << "llh " << loglh << ", beta " << (*p) * rad_to_deg << endl;
//    fval = loglh;
//}

//void MinuitFitter::DDFCN(int& npar, double* grad, double &fval,
//                         double *p, int iflag) {
//    m_pars->set_beta(*p);

//    auto pdf = Cfg::pdf();
//    DDMap csmap;  // cached coefficients
//    for (auto bbin = -8; bbin <= 8; bbin++) if (bbin)
//        for (auto dbin = -8; dbin <= 8; dbin++) if (dbin)
//            csmap[bbin][dbin] = m_pars->coefs(bbin, dbin);

//    uint64_t cnt = 0;
//    double loglh = 0.;
//    for (const auto& evt : sdevtv) {
//        if (!single_bin || (single_bin == abs(evt.Bbin()))) {
//            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
//            loglh += pdf(evt.t(), evt.Tag(), cs.first, cs.second);
//            cnt++;
//        }
//    }
//    loglh *= -2.;
//    cout << "llh: " << loglh
//         << ", beta: " << *p
////         << ", rb: " << m_pars->rb()
//         << ", cnt: " << cnt
//         << endl;
//    fval = loglh;
//}

//void MinuitFitter::ReadDBData(const string &input_file_name) {
//    sdevtv.clear();
//    ifstream ifile(input_file_name, std::ifstream::in);
//    if (!ifile.is_open()) {
//        cerr << "Can't open file " << input_file_name << endl;
//        return;
//    }
//    while (ifile.good()) {
//        static double time;
//        static int16_t tag, dbin, bbin;
//        ifile >> time >> tag >> dbin >> bbin;
//        sdevtv.emplace_back(time, tag, dbin, bbin);
//    }
//    cout << sdevtv.size() << " events found" << endl;
//    return;
//}

//void MinuitFitter::DBFit(const string &input_file_name, unique_ptr<AbsDDPars> &pars) {
//    ReadDBData(input_file_name);
//    m_pars = std::move(pars);

//    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
//    std::default_random_engine rndmeng;
//    rndmeng.seed(std::random_device {}());

//    VFit::SetDefaultFitter("Minuit");
//    const int npar = 1;
//    VFit* fitter = VFit::Fitter(0, npar);
//    fitter->SetParameter(0, "beta", udist(rndmeng), 1., -10., 10.);
//    fitter->SetFCN(DDFCN);

//    double arglist[100];
//    // Set pring level
//    arglist[0] = 0;
//    fitter->ExecuteCommand("SET PRINT", arglist, 2);
//    // Minimize
//    arglist[0] = 0;  // number of calls
//    arglist[1] = 1.e-10;  // tolerance
//    fitter->ExecuteCommand("MIGRAD", arglist, 2);

//    // Get results
//    auto beta_fit = fitter->GetParameter(0) * rad_to_deg;
//    while (beta_fit < 0) beta_fit += 180;
//    while (beta_fit > 180) beta_fit -= 180;
//    auto beta_error = fitter->GetParError(0) * rad_to_deg;
//    cout << "beta = " << beta_fit << " +- " << beta_error << endl;
//}

//void MinuitFitter::MakeFit() {
//    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
//    std::default_random_engine rndmeng;
//    rndmeng.seed(std::random_device {}());

//    VFit::SetDefaultFitter("Minuit");
//    const int npar = sfmode == fitmode::Full ? 3 : 1;
//    VFit* fitter = VFit::Fitter(0, npar);
//    fitter->SetParameter(0, "beta", udist(rndmeng), 1., 0., 0.);
//    if (sfmode == fitmode::Full) {
//        fitter->SetParameter(1, "rb", 0.02, 0.01, 0., .2);
//        fitter->SetParameter(2, "db", udist(rndmeng), 1., -2.*M_PI, 2.*M_PI);
//    }
//    fitter->SetFCN(TheFCN);

//    double arglist[100];
//    // Set pring level
//    arglist[0] = 0;
//    fitter->ExecuteCommand("SET PRINT", arglist, 2);
//    // Minimize
//    arglist[0] = 0;  // number of calls
//    arglist[1] = 1.e-10;  // tolerance
//    fitter->ExecuteCommand("MIGRAD", arglist, 2);

//    // Get results
//    double beta_fit;
//    beta_fit = fitter->GetParameter(0) * rad_to_deg;
//    while (beta_fit < 0) beta_fit += 180;
//    while (beta_fit > 180) beta_fit -= 180;
//    double beta_error = fitter->GetParError(0) * rad_to_deg;
//    cout << "beta = " << beta_fit << " +- " << beta_error << endl;

//    double rb_fit = fitter->GetParameter(1);
//    double rb_error = fitter->GetParError(1);
//    cout << "  rb = " << rb_fit << " +- " << rb_error << endl;

//    double db_fit = fitter->GetParameter(2) * rad_to_deg;
//    while (db_fit < 0) beta_fit += 360;
//    while (db_fit > 360) beta_fit -= 360;
//    double db_error = fitter->GetParError(2) * rad_to_deg;
//    cout << "delb = " << db_fit << " +- " << db_error << endl;
//}

//void MinuitFitter::ReadData(const string &input_file_name) {
//    sevtv.clear();
//    ifstream ifile(input_file_name, std::ifstream::in);
//    if (!ifile.is_open()) {
//        cerr << "Can't open file " << input_file_name << endl;
//        return;
//    }
//    double time;
//    int16_t tag;
//    int type;
//    int16_t bin;
//    while (ifile.good()) {
//        ifile >> time >> tag >> bin >> type;
//        sevtv.push_back(BEvt(time, tag, bin, static_cast<dtypes>(type)));
//    }
//    ifile.close();
//}
