#include "cpvminimizer.h"

#include <iostream>
#include <map>
#include <cmath>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mylibs/libTatami/toypdf.h"

#include "rndmd.h"
#include "cfg.h"

#include "absddpars.h"
#include "bevt.h"
#include "dbevt.h"

using std::cout;
using std::endl;
using std::string;
using std::to_string;
using std::abs;
using std::vector;
using std::unique_ptr;

using Functor = ROOT::Math::Functor;
using Factory = ROOT::Math::Factory;

using Pdf = libTatami::ToyPdf;

const std::map<fitmode, bool> CPVMinimizer::fitphase {
    {fitmode::Full,        true},
    {fitmode::Simple,      false},
    {fitmode::Corrected,   true},
    {fitmode::Approx,      true},
    {fitmode::Dh,          false},
    {fitmode::DhCorrected, false},
    {fitmode::DKs,         false}
};

// Static variables //
string CPVMinimizer::m_min_name("Minuit2");
string CPVMinimizer::m_min_alg("Migrad");
uint16_t CPVMinimizer::m_single_bin(0);
bool CPVMinimizer::m_rndm_init(false);
fitmode CPVMinimizer::m_mode(fitmode::Full);
Pdf CPVMinimizer::m_pdf(Cfg::pdf());
bool CPVMinimizer::m_dh_flag(false);
unique_ptr<AbsDDPars> CPVMinimizer::m_pars;
vector<BEvt> CPVMinimizer::sevtv;
vector<DBEvt> CPVMinimizer::sdevtv;

double CPVMinimizer::DDLH() {
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
        if (!m_single_bin || (m_single_bin == abs(evt.Bbin()))) {
            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
            auto pdfv = m_pdf(evt.t(), evt.Tag(), cs.first, cs.second);
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
        cout << "DDLH: bad/good counts " << badcnt << "/" << cnt << endl;
    return loglh;
}

double CPVMinimizer::CPLH() {
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
        if (!m_single_bin || (m_single_bin == abs(evt.Bin()))) {
            auto& cs = csmap[evt.Bin()][evt.Type()];
            auto pdfv = m_pdf(evt.t(), evt.Tag(), cs.first, cs.second);
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
        cout << "CPLH: bad/good counts " << badcnt << "/" << cnt << endl;
    return loglh;
}

double CPVMinimizer::DhLH() {
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
        if (!m_single_bin || (m_single_bin == abs(evt.Bin()))) {
            auto& cs = csmap[evt.Bin()][evt.Type()];
            auto pdfv = m_pdf(evt.t(), evt.Tag(), cs.first, cs.second);
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
        cout << "CPLH: bad/good counts " << badcnt << "/" << cnt << endl;
    return loglh;
}

double CPVMinimizer::lhfcn(const double* x) {
    m_pars->setParam("beta", *x);
    if (m_mode == fitmode::DKs)
        m_pars->setParam("gamma", *(x+1));
    if (!m_dh_flag || (m_mode != fitmode::Simple)) {
        vector<double> cvec(m_nbins);
        vector<double> svec(m_nbins);
        for (auto idx = 1; idx <= m_nbins; idx++) {
            cvec[idx] = *(x + 2*idx - 1);
            svec[idx] = *(x + 2*idx);
        }
        m_pars->set_c(cvec);
        m_pars->set_s(svec);
    }
    double fval = m_dh_flag ? -2. * DhLH() : -2. * (DDLH() + CPLH());
    cout << "llh: " << fval << ", beta: " << *x * rad_to_deg;
    if (m_mode == fitmode::DKs)
        cout << ", gamma: " << *(x+1) * rad_to_deg;
    cout << endl;
    return fval;
}

void CPVMinimizer::MakeFit(std::unique_ptr<AbsDDPars>& pars,
                           fitmode mode, uint16_t sb) {
    m_pars = std::move(pars);
    m_single_bin = sb;
    m_mode = mode;
    m_dh_flag = m_mode == fitmode::Dh ||
                m_mode == fitmode::DhCorrected ||
                m_mode == fitmode::DKs;
    RndmD rndm(-M_PI, M_PI);
    RndmD urndm(-1., 1.);
    // create minimizer giving a name and a name (optionally) for the specific
    // algorithm
    // possible choices are:
    //     minName                  algoName
    // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
    //  Minuit2                     Fumili2
    //  Fumili
    //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
    //                              BFGS2, SteepestDescent
    //  GSLMultiFit
    //   GSLSimAn
    //   Genetic
    auto* minimum = Factory::CreateMinimizer(m_min_name, m_min_alg);
    // set tolerance , etc...
    minimum->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
//    minimum->SetMaxIterations(10000);  // for GSL
    minimum->SetTolerance(0.001);
    minimum->SetPrintLevel(2);

    // create function wrapper for minimizer
    // a IMultiGenFunction type
    Functor f(&lhfcn, fitphase.at(mode) ? 18 : 2);
    minimum->SetFunction(f);

    // Set the free variables to be minimized !
    auto ibeta = m_rndm_init ? rndm() : m_pars->getParam("beta");
    auto igamm = m_rndm_init ? rndm() : m_pars->getParam("gamma");
    minimum->SetVariable(0, "beta", ibeta, 1.);
    minimum->SetVariable(1, "gamma", igamm, 1.);
    if (fitphase.at(mode)) {
        auto cvec = m_pars->get_c();
        auto svec = m_pars->get_s();
        for (auto idx = 1; idx <= m_nbins; idx++) {
            static double cini, sini;
            if (m_rndm_init) {
                cini = urndm();
                sini = urndm();
            } else {
                cini = cvec[idx-1];
                sini = svec[idx-1];
            }
            minimum->SetVariable(2 * idx, "c" + to_string(idx), cini, 0.1);
            minimum->SetVariable(2 * idx + 1, "s" + to_string(idx), sini, 0.1);
        }
    }

    // do the minimization
    minimum->Minimize();
    const double *xs = minimum->X();
    const double *xe = minimum->Errors();
    cout << "Minimum: f(" << xs[0] << "): " << minimum->MinValue()  << endl;

    cout << "beta  = " << xs[0] * rad_to_deg << " +- "
                       << xe[0] * rad_to_deg << endl;
    if (m_mode == fitmode::DKs)
         cout << "gamma = " << xs[1] * rad_to_deg << " +- "
                            << xe[1] * rad_to_deg << endl;
    if (fitphase.at(mode))
        for (auto idx = 1; idx <= m_nbins; idx++)
            cout << "c" + to_string(idx) + " = "
                 << xs[2 * idx] << " +- " << xe[2 * idx] << endl
                 << "s" + to_string(idx) + " = "
                 << xs[2 * idx + 1] << " +- " << xe[2 * idx + 1] << endl;
}

void CPVMinimizer::ReadData(const string &input_file_name) {
    if (sevtv.empty()) {
        sevtv = ReadEvents<BEvt>(input_file_name);
    } else {
        auto newev = ReadEvents<BEvt>(input_file_name);
        sevtv.insert(sevtv.end(), newev.begin(), newev.end());
    }
}

void CPVMinimizer::ReadBData(const string &input_file_name) {
    if (sevtv.empty()) {
        sdevtv = ReadEvents<DBEvt>(input_file_name);
    } else {
        auto newev = ReadEvents<DBEvt>(input_file_name);
        sdevtv.insert(sdevtv.end(), newev.begin(), newev.end());
    }
}

void CPVMinimizer::FlushData() {
    sevtv.clear();
    sdevtv.clear();
}
