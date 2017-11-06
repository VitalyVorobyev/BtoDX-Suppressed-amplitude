//#ifndef MINUITFITTER_H
//#define MINUITFITTER_H

//#include <vector>
//#include <string>
//#include <cstdint>
//#include <memory>
//#include <fstream>
//#include <random>  // std::default_random_engine
//#include <cmath>
//#include <map>
//#include <iostream>

//#include "TVirtualFitter.h"

//#include "dbevt.h"
//#include "fitmodes.h"
//#include "absddpars.h"
//#include "cfg.h"
//#include "lambda.h"
//#include "binnedparams.h"

//template<class P>
//class MinuitFitter {
//    static constexpr double rad_to_deg = 180. / M_PI;
//    static constexpr int m_nbins = 8;
//    /** B0 -> Dcp pi+ pi- events */
//    static std::vector<std::unique_ptr<BEvt>> sevtv;
//    /** B0 -> (Ks0 pi+ pi-)_D0 pi+ pi- events */
//    static std::vector<std::unique_ptr<DBEvt>> sdevtv;
//    static uint16_t single_bin;
//    static fitmode sfmode;
//    static void TheFCN(int& npar, double* grad,
//                       double &fval, double *p, int iflag);
//    static void DDFCN(int& npar, double* grad,
//                      double &fval, double *p, int iflag);
//    static void FullDDFCN(int& npar, double* grad,
//                          double &fval, double *p, int iflag);
//    static std::unique_ptr<P> m_pars;

//public:
//    explicit MinuitFitter(fitmode mode = fitmode::Simple) {
//        sfmode = mode;
//    }
//    MinuitFitter(const std::string& input_file_name,
//                 fitmode mode = fitmode::Simple) : MinuitFitter(mode) {
//        ReadData(input_file_name);
//    }
//    MinuitFitter(const std::vector<std::string>& input_files,
//                 fitmode mode = fitmode::Simple) : MinuitFitter(mode) {
//        for (const auto& fname : input_files) ReadData(fname);
//    }
//    /**
//     * @brief ReadData. Read event vector from text file
//     * @param input_file_name. Input text file
//     */
//    void ReadData(const std::string &input_file_name);
//    /**
//     * @brief TheFCN. Unbinned Log Likelihood Function
//     * @param npar number of parameters to fit
//     * @param grad precomputed gradient (optional)
//     * @param fval fcn value (to be filled)
//     * @param p parameters vector
//     * @param iflag flag
//     */
//    void MakeFit();
//    void ReadDBData(const std::string &input_file_name);
//    void FullDBFit(std::unique_ptr<P> &pars);
//    void DBFit(const std::string &input_file_name,
//               std::unique_ptr<P>& pars);
//    void setSingleBin(uint16_t sb) {single_bin = sb;}
//    void unsetSingleBin() {setSingleBin(0);}
//};

//template<class P>
//std::vector<std::unique_ptr<BEvt>> MinuitFitter<P>::sevtv;
//template<class P>
//std::vector<std::unique_ptr<DBEvt>> MinuitFitter<P>::sdevtv;
//template<class P>
//uint16_t MinuitFitter<P>::single_bin;
//template<class P>
//fitmode MinuitFitter<P>::sfmode;
//template<class P>
//std::unique_ptr<P> MinuitFitter<P>::m_pars;

//template<class P>
//void MinuitFitter<P>::TheFCN(int& npar, double* grad,
//                          double &fval, double *p, int iflag) {
//    static auto spdf = Cfg::pdf();
//    static auto bpars = Cfg::bpars(sfmode == fitmode::Approx);
//    static std::unordered_map<dtypes, Lambda_b2cud> slambdas;
//    if (slambdas.empty()) {  // initialize once
//        if (sfmode != fitmode::Simple) {
//            slambdas.emplace(std::make_pair(dtypes::CPp, Cfg::lamf(dtypes::CPp)));
//            slambdas.emplace(std::make_pair(dtypes::CPn, Cfg::lamf(dtypes::CPn)));
//            slambdas.emplace(std::make_pair(dtypes::KPI, Cfg::lamf(dtypes::KPI)));
//            slambdas.emplace(std::make_pair(dtypes::PIK, Cfg::lamf(dtypes::PIK)));
//        } else {
//            bpars->setRb(0.);
//            slambdas.emplace(std::make_pair(dtypes::CPp, Cfg::lamf0(dtypes::CPp)));
//            slambdas.emplace(std::make_pair(dtypes::CPn, Cfg::lamf0(dtypes::CPn)));
//            slambdas.emplace(std::make_pair(dtypes::KPI, Cfg::lamf0(dtypes::KPI)));
//            slambdas.emplace(std::make_pair(dtypes::PIK, Cfg::lamf0(dtypes::PIK)));
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
//        if (single_bin && (abs(evt->Bin()) != single_bin)) continue;
//        double ccoef, scoef;
//        if (evt->Type() == dtypes::KsPIPI) {
//            auto coefs = bpars->coefs(evt->Bin());
//            ccoef = coefs.first;
//            scoef = coefs.second;
//        } else {
//            auto &lamf = slambdas.find(evt->Type())->second;
//            ccoef = lamf.ccoef();
//            scoef = lamf.scoef();
//        }
//        loglh += -2. * log(spdf(evt->t(), evt->Tag(), ccoef, scoef));
//    }
////    std::cout << "llh " << loglh << ", beta " << (*p) * rad_to_deg << std::endl;
//    fval = loglh;
//}

//template<class P>
//void MinuitFitter<P>::DDFCN(int& npar, double* grad, double &fval,
//                            double *p, int iflag) {
//    auto pdf = Cfg::pdf();
//    std::unordered_map<int16_t,
//            std::unordered_map<int16_t,
//            std::pair<double, double>>> csmap;  // cached coefficients
//    m_pars->set_beta(*p);
//    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin)
//        for (auto dbin = -m_nbins; dbin <= m_nbins; dbin++) if (dbin)
//            csmap[bbin][dbin] = m_pars->coefs(bbin, dbin);

//    uint64_t cnt = 0;
//    uint64_t badcnt = 0;
//    double loglh = 0.;
//    for (const auto& evt : sdevtv) {
//        if (!single_bin || (single_bin == abs(evt.Bbin()))) {
//            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
//            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
//            if (pdfv > 0) {
//                loglh += std::log(pdfv);
//                cnt++;
//            } else {
//                loglh -= 100. * pdfv;
//                badcnt++;
//            }
//        }
//    }
//    loglh *= -2.;
//    std::cout << "llh: " << loglh
//         << ", beta: " << *p
//         << ", cnt: " << cnt
//         << ", badcnt: " << badcnt
//         << std::endl;
//    fval = loglh;
//}

//template<class P>
//double MinuitFitter<P>::DDLH(const libTatami::ToyPdf& pdf) {
//    if (sdevtv.empty()) return 0.;
//    // Cache coefficients
//    MapCS csmap;  // cached coefficients
//    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin)
//        for (auto dbin = -m_nbins; dbin <= m_nbins; dbin++) if (dbin)
//            csmap[bbin][dbin] = m_pars->coefs(bbin, dbin, dtypes::KsPIPI);
//    // Calculate loglh
//    uint64_t cnt = 0;
//    uint64_t badcnt = 0;
//    double loglh = 0.;
//    for (const auto& evt : sdevtv) {
//        if (!single_bin || (single_bin == abs(evt.Bbin()))) {
//            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
//            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
//            if (pdfv > 0) {
//                loglh += std::log(pdfv);
//                cnt++;
//            } else {
//                loglh -= 100. * pdfv;
//                badcnt++;
//            }
//        }
//    }
//    if (badcnt)
//        cout << "DDLH: bad/good counts " << badcnt << "/" << cnt << endl;
//    return loglh;
//}

//template<class P>
//double MinuitFitter<P>::CPLH(const libTatami::ToyPdf& pdf) {
//    if (sevtv.empty()) return 0.;
//    SMapCS scsmap;  // cached coefficients for B0 -> Dcp pi+ pi-
//    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin) {
//        scsmap[bbin][dtypes::CPp] = m_pars->coefs(bbin, 0, dtypes::CPp);
//        scsmap[bbin][dtypes::CPn] = m_pars->coefs(bbin, 0, dtypes::CPn);
//    }
//    // Calculate loglh
//    uint64_t cnt = 0;
//    uint64_t badcnt = 0;
//    double loglh = 0.;
//    for (const auto& evt : sevtv) {

//    }
//    if (badcnt)
//        cout << "CPLH: bad/good counts " << badcnt << "/" << cnt << endl;
//    return loglh;
//}

//template<class P>
//void MinuitFitter<P>::FullFCN(int& npar, double* grad, double &fval,
//                         double *p, int iflag) {
//    // Update parameters
//    std::vector<double> cvec(m_nbins);
//    std::vector<double> svec(m_nbins);
//    for (auto idx = 0; idx < m_nbins; idx++) {
//        cvec[idx] = *(p + 2*idx + 1);
//        svec[idx] = *(p + 2*idx + 2);
//    }
//    m_pars->set_beta(*p);
//    m_pars->set_c(cvec);
//    m_pars->set_s(svec);

//    auto pdf = Cfg::pdf();
//    const auto ddlh = MinuitFitter<P>::DDLH(pdf);
//    const auto cplh = MinuitFitter<P>::CPLH(pdf);

//    fval = -2. * (ddlh + cplh);
//    std::cout << "llh: " << loglh << ", beta: " << *p << ", cnt: " << cnt
//         << ", badcnt: " << badcnt << std::endl;
//    fval = loglh;
//}

//template<class P>
//void MinuitFitter<P>::FullDBFit(std::unique_ptr<P> &pars) {
//    m_pars = std::move(pars);

//    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
//    std::uniform_real_distribution<double> uunit(-1., 1.);
//    std::default_random_engine rndmeng;
//    rndmeng.seed(std::random_device {}());

//    TVirtualFitter::SetDefaultFitter("Minuit");
//    const int npar = 1;
//    TVirtualFitter* fitter = TVirtualFitter::Fitter(0, npar);
//    fitter->SetParameter(0, "beta", udist(rndmeng), 1., -10., 10.);
//    for (auto idx = 1; idx <= m_nbins; idx++) {
//        static std::string name;
//        name = "c" + std::to_string(idx);
//        fitter->SetParameter(2 * idx - 1, name.c_str(), uunit(rndmeng), 0.5, -1.1, 1.1);
//        name = "s" + std::to_string(idx);
//        fitter->SetParameter(2 * idx, name.c_str(), uunit(rndmeng), 0.5, -1.1, 1.1);
//    }
//    fitter->SetFCN(FullDDFCN);

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
//    std::cout << "beta = " << beta_fit << " +- " << beta_error << std::endl;

//    for (auto idx = 0; idx < m_nbins; idx++) {
//        std::cout << "c" + std::to_string(idx+1) + " = "
//                  << fitter->GetParameter(2 * idx + 1) << " +- "
//                  << fitter->GetParError(2 * idx + 1) << std::endl;
//        std::cout << "s" + std::to_string(idx+1) + " = "
//                  << fitter->GetParameter(2 * idx + 2) << " +- "
//                  << fitter->GetParError(2 * idx + 2) << std::endl;
//    }
//}

//template<class P>
//void MinuitFitter<P>::DBFit(const std::string &input_file_name,
//                         std::unique_ptr<P> &pars) {
//    ReadDBData(input_file_name);
//    m_pars = std::move(pars);

//    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
//    std::default_random_engine rndmeng;
//    rndmeng.seed(std::random_device {}());

//    TVirtualFitter::SetDefaultFitter("Minuit");
//    const int npar = 1;
//    TVirtualFitter* fitter = TVirtualFitter::Fitter(0, npar);
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
//    std::cout << "beta = " << beta_fit << " +- " << beta_error << std::endl;
//}

//template<class P>
//void MinuitFitter<P>::MakeFit() {
//    std::uniform_real_distribution<double> udist(-M_PI, M_PI);
//    std::default_random_engine rndmeng;
//    rndmeng.seed(std::random_device {}());

//    TVirtualFitter::SetDefaultFitter("Minuit");
//    const int npar = sfmode == fitmode::Full ? 3 : 1;
//    TVirtualFitter* fitter = TVirtualFitter::Fitter(0, npar);
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
//    std::cout << "beta = " << beta_fit << " +- " << beta_error << std::endl;

//    double rb_fit = fitter->GetParameter(1);
//    double rb_error = fitter->GetParError(1);
//    std::cout << "  rb = " << rb_fit << " +- " << rb_error << std::endl;

//    double db_fit = fitter->GetParameter(2) * rad_to_deg;
//    while (db_fit < 0) beta_fit += 360;
//    while (db_fit > 360) beta_fit -= 360;
//    double db_error = fitter->GetParError(2) * rad_to_deg;
//    std::cout << "delb = " << db_fit << " +- " << db_error << std::endl;
//}

//template<class P>
//void MinuitFitter<P>::ReadData(const std::string &input_file_name) {
//    sevtv.clear();
//    std::ifstream ifile(input_file_name, std::ifstream::in);
//    if (!ifile.is_open()) {
//        std::cerr << "Can't open file " << input_file_name << std::endl;
//        return;
//    }
//    cout << "MinuitFitter: read data from " << input_file_name << endl;
//    double time;
//    int16_t tag;
//    int type;
//    int16_t bin;
//    while (ifile.good()) {
//        ifile >> time >> tag >> bin >> type;
//        sevtv.push_back(BEvt(time, tag, bin, static_cast<dtypes>(type)));
//    }
//    cout << sevtv.size() << " events found" << endl;
//    ifile.close();
//}

//template<class P>
//void MinuitFitter<P>::ReadDBData(const std::string &input_file_name) {
//    sdevtv.clear();
//    std::ifstream ifile(input_file_name, std::ifstream::in);
//    if (!ifile.is_open()) {
//        std::cerr << "Can't open file " << input_file_name << std::endl;
//        return;
//    }
//    cout << "MinuitFitter: read data from " << input_file_name << endl;
//    while (ifile.good()) {
//        static double time;
//        static int16_t tag, dbin, bbin;
//        ifile >> time >> tag >> dbin >> bbin;
//        sdevtv.emplace_back(time, tag, dbin, bbin);
//    }
//    std::cout << sdevtv.size() << " events found" << std::endl;
//    return;
//}

//#endif // MINUITFITTER_H
