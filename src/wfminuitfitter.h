#ifndef WFMINUITFITTER_H
#define WFMINUITFITTER_H

#include <vector>
#include <string>
#include <cstdint>
#include <memory>
#include <fstream>
#include <cmath>
#include <map>
#include <iostream>
#include <exception>

#include "TVirtualFitter.h"

#include "bevt.h"
#include "dbevt.h"
#include "fitmodes.h"
#include "absddpars.h"
#include "cfg.h"
#include "binnedparams.h"
#include "rndmd.h"

template<class P>
class WFMinuitFitter {
    static const std::map<fitmode, bool> fitphase;
    /** Binnned double Dalitz C and S keeper */
    using MapCS = std::unordered_map<int16_t,
                  std::unordered_map<int16_t,
                  std::pair<double, double>>>;
    /** Binnned Dalitz C and S keeper */
    using SMapCS = std::unordered_map<int16_t,
                   std::unordered_map<dtypes,
                   std::pair<double, double>>>;
    static constexpr double rad_to_deg = 180. / M_PI;
    static constexpr int m_nbins = 8;
    static fitmode m_mode;
    /** B0 -> Dcp pi+ pi- events */
    static std::vector<BEvt> sevtv;
    /** B0 -> (Ks0 pi+ pi-)_D0 pi+ pi- events */
    static std::vector<DBEvt> sdevtv;
    static uint16_t single_bin;
    static bool m_rndm_init;
    /** @brief Make fit with only angle beta as a free parameter */
    static void SimpleFCN(int& npar, double* grad,
                          double &fval, double *p, int iflag);
    /** @brief Complete fit procesure with andgle beta and
     *  16 phase parametres free */
    static void FullFCN(int& npar, double* grad,
                        double &fval, double *p, int iflag);
    /** P is DDBPars or DDMPars */
    static std::unique_ptr<P> m_pars;
    static double DDLH(libTatami::ToyPdf& pdf);
    static double CPLH(libTatami::ToyPdf &pdf);
    static double DhLH(libTatami::ToyPdf& pdf);
    /** Private constructor */
    WFMinuitFitter() {}

 public:
    static void ReadData(const std::string &input_file_name) {
        if (sevtv.empty()) {
            sevtv = ReadEvents<BEvt>(input_file_name);
        } else {
            auto newev = ReadEvents<BEvt>(input_file_name);
            sevtv.insert(sevtv.end(), newev.begin(), newev.end());
        }
    }
    static void ReadBData(const std::string &input_file_name) {
        if (sevtv.empty()) {
            sdevtv = ReadEvents<DBEvt>(input_file_name);
        } else {
            auto newev = ReadEvents<DBEvt>(input_file_name);
            sdevtv.insert(sdevtv.end(), newev.begin(), newev.end());
        }
    }
    static void FlushData() {
        sevtv.clear();
        sdevtv.clear();
    }
    static void TheFCN(int& npar, double* grad, double &fval,
                       double *p, int iflag);
    static void TheFCNDh(int& npar, double* grad, double &fval,
                         double *p, int iflag);
    static void MakeFit(std::unique_ptr<P>& pars, fitmode mode, uint16_t sb=0);
    static void RndmInit(bool x=true) {m_rndm_init = x;}
};

template<class P> std::vector<DBEvt> WFMinuitFitter<P>::sdevtv;
template<class P> std::vector<BEvt> WFMinuitFitter<P>::sevtv;
template<class P> uint16_t WFMinuitFitter<P>::single_bin;
template<class P> fitmode WFMinuitFitter<P>::m_mode;
template<class P> std::unique_ptr<P> WFMinuitFitter<P>::m_pars;
template<class P> bool WFMinuitFitter<P>::m_rndm_init = false;
template<class P> const std::map<fitmode, bool> WFMinuitFitter<P>::fitphase {
    {fitmode::Full,        true},
    {fitmode::Simple,      false},
    {fitmode::Corrected,   true},
    {fitmode::Approx,      true},
    {fitmode::Dh,          false},
    {fitmode::DhCorrected, false}
};

template<class P>
double WFMinuitFitter<P>::DDLH(libTatami::ToyPdf& pdf) {
    if (sdevtv.empty()) return 0.;
    // Cache coefficients
    MapCS csmap;  // cached coefficients
    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin)
        for (auto dbin = -m_nbins; dbin <= m_nbins; dbin++) if (dbin)
            csmap[bbin][dbin] = m_pars->coefs(bbin, dbin, dtypes::KsPIPI);
    // Calculate loglh
    uint64_t cnt = 0;
    uint64_t badcnt = 0;
    double loglh = 0.;
    for (const auto& evt : sdevtv) {
        if (!single_bin || (single_bin == abs(evt.Bbin()))) {
            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv))
                throw std::runtime_error("nan pdf");
            if (pdfv > 0) {
                loglh += std::log(pdfv);
                cnt++;
            } else {
                loglh -= 1000. * pdfv;
                badcnt++;
            }
        }
    }
    if (badcnt || !cnt)
        std::cout << "DDLH: bad/good counts " << badcnt << "/"
                  << cnt << std::endl;
    return loglh;
}

template<class P>
double WFMinuitFitter<P>::CPLH(libTatami::ToyPdf& pdf) {
    if (sevtv.empty()) return 0.;
    SMapCS csmap;  // cached coefficients for B0 -> Dcp pi+ pi-
    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin) {
        csmap[bbin][dtypes::CPp] = m_pars->coefs(bbin, 0, dtypes::CPp);
        csmap[bbin][dtypes::CPn] = m_pars->coefs(bbin, 0, dtypes::CPn);
    }
    // Calculate loglh
    uint64_t cnt = 0;
    uint64_t badcnt = 0;
    double loglh = 0.;
    for (const auto& evt : sevtv) {
        if (!single_bin || (single_bin == abs(evt.Bin()))) {
            auto& cs = csmap[evt.Bin()][evt.Type()];
            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv))
                throw std::runtime_error("nan pdf");
            if (pdfv > 0) {
                loglh += std::log(pdfv);
                cnt++;
            } else {
                loglh -= 1000. * pdfv;
                badcnt++;
            }
        }
    }
    if (badcnt || !cnt)
        std::cout << "CPLH: bad/good counts " << badcnt
                  << "/" << cnt << std::endl;
    return loglh;
}

template<class P>
double WFMinuitFitter<P>::DhLH(libTatami::ToyPdf& pdf) {
    if (sevtv.empty()) return 0.;
    SMapCS csmap;  // cached coefficients for B0 -> Dcp h0
    csmap[0][dtypes::DhCPp] = m_pars->coefs(0, 0, dtypes::DhCPp);
    csmap[0][dtypes::DhCPn] = m_pars->coefs(0, 0, dtypes::DhCPn);
    for (auto dbin = -m_nbins; dbin <= m_nbins; dbin++) if (dbin)
        csmap[dbin][dtypes::Dh] = m_pars->coefs(0, dbin, dtypes::Dh);
    // Calculate loglh
    uint64_t cnt = 0;
    uint64_t badcnt = 0;
    double loglh = 0.;
    for (const auto& evt : sevtv) {
        if (!single_bin || (single_bin == abs(evt.Bin()))) {
            auto& cs = csmap[evt.Bin()][evt.Type()];
            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv))
                throw std::runtime_error("nan pdf");
            if (pdfv > 0) {
                loglh += std::log(pdfv);
                cnt++;
            } else {
                loglh -= 1000. * pdfv;
                badcnt++;
            }
        }
    }
    if (badcnt || !cnt)
        std::cout << "CPLH: bad/good counts " << badcnt
                  << "/" << cnt << std::endl;
    return loglh;
}

template<class P>
void WFMinuitFitter<P>::TheFCN(int& npar, double* grad, double &fval,
                              double *p, int iflag) {
    // Update parameters
    m_pars->set_beta(*p);
    if (m_mode != fitmode::Simple) {
        std::vector<double> cvec(m_nbins);
        std::vector<double> svec(m_nbins);
        for (auto idx = 0; idx < m_nbins; idx++) {
            cvec[idx] = *(p + 2*idx + 1);
            svec[idx] = *(p + 2*idx + 2);
        }
        m_pars->set_c(cvec);
        m_pars->set_s(svec);
    }
    auto pdf = Cfg::pdf();
    const auto ddlh = WFMinuitFitter<P>::DDLH(pdf);
    const auto cplh = WFMinuitFitter<P>::CPLH(pdf);
    fval = -2. * (ddlh + cplh);
    std::cout << "llh: " << fval << ", beta: " << *p << std::endl;
}

template<class P>
void WFMinuitFitter<P>::TheFCNDh(int& npar, double* grad, double &fval,
                                 double *p, int iflag) {
    // Update parameters
    m_pars->set_beta(*p);
    auto pdf = Cfg::pdf();
    fval = -2. * WFMinuitFitter<P>::DhLH(pdf);
    std::cout << "llh: " << fval << ", beta: " << *p << std::endl;
}

template<class P>
void WFMinuitFitter<P>::MakeFit(std::unique_ptr<P> &pars, fitmode mode, uint16_t sb) {
    m_pars = std::move(pars);
    single_bin = sb;
    m_mode = mode;
    RndmD rndm(-M_PI, M_PI);
    RndmD urndm(-1., 1.);

    TVirtualFitter::SetDefaultFitter("Minuit");
    const int npar = 1;
    TVirtualFitter* fitter = TVirtualFitter::Fitter(0, npar);
    auto ibeta = m_rndm_init ? rndm() : m_pars->beta();
    fitter->SetParameter(0, "beta", ibeta, 1., -10., 10.);
    if (fitphase.at(mode)) {
        auto cvec = m_pars->get_c();
        auto svec = m_pars->get_s();
        for (auto idx = 1; idx <= m_nbins; idx++) {
            static std::string name;
            static double cini, sini;
            name = "c" + std::to_string(idx);
            if (m_rndm_init) {
                cini = urndm();
                sini = urndm();
            } else {
                cini = cvec[idx-1];
                sini = svec[idx-1];
            }
            fitter->SetParameter(2 * idx - 1, name.c_str(), cini, 0.1, -1.1, 1.1);
            name = "s" + std::to_string(idx);
            fitter->SetParameter(2 * idx, name.c_str(), sini, 0.1, -1.1, 1.1);
        }
    }
    if (mode == fitmode::Dh || mode == fitmode::DhCorrected)
        fitter->SetFCN(TheFCNDh);
    else
        fitter->SetFCN(TheFCN);
    double arglist[100];
    // Set pring level
    arglist[0] = 0;
    fitter->ExecuteCommand("SET PRINT", arglist, 2);
    // Minimize
    arglist[0] = 0;  // number of calls
    arglist[1] = 1.e-10;  // tolerance
    fitter->ExecuteCommand("MIGRAD", arglist, 2);

    // Get results
    auto beta_fit = fitter->GetParameter(0) * rad_to_deg;
    while (beta_fit < 0) beta_fit += 180;
    while (beta_fit > 180) beta_fit -= 180;
    auto beta_error = fitter->GetParError(0) * rad_to_deg;
    std::cout << "beta = " << beta_fit << " +- " << beta_error << std::endl;
    if (mode != fitmode::Simple) {
        for (auto idx = 0; idx < m_nbins; idx++) {
            std::cout << "c" + std::to_string(idx+1) + " = "
                      << fitter->GetParameter(2 * idx + 1) << " +- "
                      << fitter->GetParError(2 * idx + 1) << std::endl;
            std::cout << "s" + std::to_string(idx+1) + " = "
                      << fitter->GetParameter(2 * idx + 2) << " +- "
                      << fitter->GetParError(2 * idx + 2) << std::endl;
        }
    }
}

#endif // WFMINUITFITTER_H
