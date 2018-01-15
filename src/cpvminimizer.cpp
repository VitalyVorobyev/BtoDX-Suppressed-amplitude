#include "cpvminimizer.h"

#include <iostream>
#include <map>
#include <cmath>
#include <exception>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include "mylibs/libTatami/toypdf.h"

#include "rndmd.h"

#include "absddpars.h"
#include "bevt.h"
#include "dbevt.h"

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::to_string;
using std::abs;
using std::vector;
using std::unique_ptr;

using Functor = ROOT::Math::Functor;
using Factory = ROOT::Math::Factory;
using Fitter = ROOT::Math::Minimizer;

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
std::unique_ptr<Pdf> CPVMinimizer::m_pdf;
bool CPVMinimizer::m_dh_flag(false);
unique_ptr<AbsDDPars> CPVMinimizer::m_pars;
vector<BEvt> CPVMinimizer::sevtv;
vector<DBEvt> CPVMinimizer::sdevtv;
CPVMinimizer::ParsLookUp CPVMinimizer::m_parsLookUp;
std::set<std::string> CPVMinimizer::m_fixed_vars;
uint16_t CPVMinimizer::m_parIdx;
Cfg::ExpSetup CPVMinimizer::m_expCfg = Cfg::ExpSetup::Perfect;

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
            auto pdfv = (*m_pdf)(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv)) {
                cerr << "nan pdf" << endl;
                throw std::runtime_error("nan pdf");
            }
            if (pdfv > 0) {
                loglh += std::log(pdfv);
                cnt++;
            } else {
                loglh = -10. * pdfv;
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
            auto pdfv = (*m_pdf)(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv)) {
                cerr << "nan pdf" << endl;
                throw std::runtime_error("nan pdf");
            }
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
            auto pdfv = (*m_pdf)(evt.t(), evt.Tag(), cs.first, cs.second);
            if (std::isnan(pdfv)) {
                cerr << "nan pdf" << endl;
                throw std::runtime_error("nan pdf");
            }
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
    if (m_fixed_vars.find("beta") == m_fixed_vars.end())
        m_pars->setParam("beta", *x);
    if (m_fixed_vars.find("gamma") == m_fixed_vars.end())
        m_pars->setParam("gamma", *(x+1));
    if (m_fixed_vars.find("rb") == m_fixed_vars.end())
        m_pars->setParam("rb", *(x+2));
    if (m_fixed_vars.find("delb") == m_fixed_vars.end())
        m_pars->setParam("delb", *(x+3));
    if (!m_dh_flag || (m_mode != fitmode::Simple)) {
        vector<double> cvec(m_nbins);
        vector<double> svec(m_nbins);
        for (auto idx = 1; idx <= m_nbins; idx++) {
            cvec[idx] = *(x + 2*idx + 1);
            svec[idx] = *(x + 2*idx + 2);
        }
        m_pars->set_c(cvec);
        m_pars->set_s(svec);
    }
    double fval = m_dh_flag ? -2. * DhLH() : -2. * (DDLH() + CPLH());
    cout << "llh: " << fval << ", beta: " << *x * rad_to_deg;
    if (m_mode == fitmode::DKs)
        cout << ", gamma: " << *(x+1) * rad_to_deg
             << ", rb: " << *(x+2)
             << ", delb: " << *(x+3) * rad_to_deg
             << endl;
    return fval;
}

uint16_t CPVMinimizer::addPar(Fitter& f,
                              const std::string& name,
                              double init, double var) {
    cout << "CPVMinimizer: new parameter " << name << endl;
    f.SetVariable(m_parIdx, name, init, var);
    m_parsLookUp.emplace(name, m_parIdx++);
    if (m_parsLookUp.size() != m_parIdx) {
        cerr << "Wrong parameters counting: " << m_parsLookUp.size()
             << " vs. " << m_parIdx << endl;
        throw new std::runtime_error("Wrong parameters counting");
    }
    return m_parIdx;
}

uint16_t CPVMinimizer::fixPars(Fitter& f) {
    for (const auto& name : m_fixed_vars) {
        const auto& var = m_parsLookUp.find(name);
        if (var != m_parsLookUp.end()) {
            cout << "CPVMinimizer: fixing " << name
                 << "(" << var->second << ")" << endl;
            f.FixVariable(var->second);
        } else {
            cerr << "can't find and fix variable" << name << endl;
            throw new std::runtime_error("Wrong variable name " + name);
        }
    }
    return m_fixed_vars.size();
}

void CPVMinimizer::fixPar(const std::string& name) {
    m_fixed_vars.emplace(name);
}

void CPVMinimizer::MakeFit(std::unique_ptr<AbsDDPars>& pars,
                           fitmode mode, uint16_t sb) {
    m_pdf = Cfg::pdf(m_expCfg);
    m_pars = std::move(pars);
    m_single_bin = sb;
    m_mode = mode;
    m_dh_flag = m_mode == fitmode::Dh ||
                m_mode == fitmode::DhCorrected ||
                m_mode == fitmode::DKs;
    m_parsLookUp.clear();
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
    auto fitter = std::unique_ptr<Fitter>(
                Factory::CreateMinimizer(m_min_name, m_min_alg));
    // set tolerance , etc...
    fitter->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
    fitter->SetTolerance(0.001);
    fitter->SetPrintLevel(2);

    // create function wrapper for minimizer
    // a IMultiGenFunction type
    Functor f(&lhfcn, 4);
    fitter->SetFunction(f);

    // Set the free variables to be minimized !
    m_parIdx = 0;
    addPar(*fitter, "beta", m_rndm_init ? rndm() : m_pars->getParam("beta"), 1.);
    addPar(*fitter, "gamma", m_pars->getParam("gamma"), 1.);
    addPar(*fitter, "rb", m_pars->getParam("rb"), 0.1);
    addPar(*fitter, "delb", m_pars->getParam("delb"), 1);
    if (fitphase.at(mode)) {
        for (auto idx = 1; idx <= m_nbins; idx++) {
            double cini = m_rndm_init ? urndm() : m_pars->get_c()[idx-1];
            double sini = m_rndm_init ? urndm() : m_pars->get_s()[idx-1];
            addPar(*fitter, "c" + to_string(idx), cini, 0.1);
            addPar(*fitter, "s" + to_string(idx), sini, 0.1);
        }
    }
    // fix some parameters
    fixPars(*fitter);
    // do the minimization
    fitter->Minimize();
    const double *xs = fitter->X();
    const double *xe = fitter->Errors();
    cout << "Minimum: f(" << xs[0] << "): " << fitter->MinValue()  << endl;

    cout << "beta  = " << xs[0] * rad_to_deg << " +- "
                       << xe[0] * rad_to_deg << endl;
    cout << "gamma = " << xs[1] * rad_to_deg << " +- "
         << xe[1] * rad_to_deg << endl
         << "   rb = " << xs[2] << " +- " << xe[2] << endl
         << " delb = " << xs[3] * rad_to_deg << " +- "
         << xe[3] * rad_to_deg << endl;
    if (fitphase.at(mode))
        for (auto idx = 1; idx <= m_nbins; idx++)
            cout << "c" + to_string(idx) + " = "
                 << xs[2 * idx + 2] << " +- " << xe[2 * idx + 2] << endl
                 << "s" + to_string(idx) + " = "
                 << xs[2 * idx + 3] << " +- " << xe[2 * idx + 3] << endl;
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
