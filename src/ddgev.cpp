#include "ddgev.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>
#include <exception>

#include "mylibs/libTatami/toypdf.h"
#include "mylibs/libTatami/toypdfgen.h"

#include "dbevt.h"
#include "ddbpars.h"
#include "ddmpars.h"
#include "absddpars.h"

using std::ofstream;
using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;

using std::move;
using std::unique_ptr;
using std::make_unique;
using std::make_move_iterator;

constexpr bool dump = false;

Cfg::ExpSetup DDGev::m_expCfg = Cfg::ExpSetup::Perfect;

using DDMap = std::map<int16_t, std::map<int16_t, double>>;
using DMap = std::map<int16_t, double>;

bool isgood_cs(const pair<double, double>& cs) {
    if (std::isnan(cs.first)     || std::isnan(cs.second) ||
        std::fabs(cs.first) > 1. || std::fabs(cs.second) > 1.) {
        cout << "check_cs: ccoef: " << cs.first
             << ", scoef: " << cs.second << endl;
        return false;
    }
    return true;
}

unique_ptr<DDGev::pBEvtVec>
DDGev::CPDalitz(uint64_t nevt, const AbsDDPars& pars, dtypes type) {
    auto pdf = Cfg::pdf(m_expCfg);  // Toy PDF
    double rate_norm = 0;
    DMap rates;
    const double alpha = pdf->alpha();
    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin) {
        auto ud = pars.rawud(bbin, 0, type);
        auto rate = ud.first + alpha * ud.second;
        rate_norm += rate;
        rates[bbin] = rate;
    }
    cout << "DDGev: " << endl
         << "  alpha: " << alpha << endl
         << "  rnorm: " << rate_norm << endl;
    libTatami::ToyPdfGen gen(*pdf);  // generator
    auto evtv = make_unique<pBEvtVec>();  // vector of events
    for (auto bbin = -m_nbins; bbin <= m_nbins; bbin++) if (bbin) {
        auto cs = pars.coefs(bbin, 0, type);
        if (!isgood_cs(cs)) {
            cerr << "DDGev::CPDalitz bbin: " << bbin << endl;
            throw new std::runtime_error("Bad C or S");
        }
        pdf->SetCS(cs);
        auto half_nevt = static_cast<uint64_t>
                (0.5 * nevt * rates[bbin] / rate_norm);
        pdf->SetTag(1);
        auto times = gen.Generate(half_nevt, true);
        for (auto t : times)
            evtv->emplace_back(make_unique<BEvt>(t, 1, bbin, type));
        pdf->SetTag(-1);
        times = gen.Generate(half_nevt, true);
        for (auto t : times)
            evtv->emplace_back(make_unique<BEvt>(t, -1, bbin, type));
    }
    return evtv;
}

void DDGev::CPDalitz(uint64_t ncpp, uint64_t ncpn, const AbsDDPars& pars,
                     const string& fname) {
    auto evtv = CPDalitz(ncpp, pars, dtypes::CPp);
    auto evtcpn = CPDalitz(ncpn, pars, dtypes::CPn);
    evtv->insert(evtv->end(),
                 make_move_iterator(evtcpn->begin()),
                 make_move_iterator(evtcpn->end()));
    cout << evtv->size() << " / " << ncpp + ncpn
         << " events generated" << endl;
    // Random shuffle positive and negative tags
    std::random_shuffle(evtv->begin(), evtv->end());
    // Write to a text file
    ofstream ofile(fname, ofstream::out);
    for (const auto& evt : *evtv)
        ofile << *evt << endl;
    ofile.close();
}

unique_ptr<DDGev::pBEvtVec>
DDGev::DhCP(uint64_t nevt, const AbsDDPars& pars, dtypes type) {
    auto pdf = Cfg::pdf(m_expCfg);  // Toy PDF
    cout << "DhCP: " << endl
         << "  alpha: " << pdf->alpha() << endl;
    libTatami::ToyPdfGen gen(*pdf);  // generator
    auto evtv = make_unique<pBEvtVec>();  // vector of events
    cout << "DhCP: getting parameters..." << endl;
    auto cs = pars.coefs(0, 0, type);
    cout << "C " << cs.first << ", S " << cs.second << endl;
    cout << "DhCP: getting parameters... done" << endl;
    if (!isgood_cs(cs)) {
        cerr << "DDGev::DhCP" << endl;
        throw new std::runtime_error("Bad C or S");
    }
    pdf->SetCS(cs);
    auto half_nevt = nevt / 2;
    pdf->SetTag(1);
    cout << "DhCP: generating..." << endl;
    auto times = gen.Generate(half_nevt, true);
    cout << "DhCP: generating... done " << times.size()
         << "/" << half_nevt << " " << nevt << endl;
    for (auto t : times)
        evtv->emplace_back(make_unique<BEvt>(t, 1, 0, type));
    pdf->SetTag(-1);
    times = gen.Generate(half_nevt, true);
    for (auto t : times)
        evtv->emplace_back(make_unique<BEvt>(t, -1, 0, type));
    return std::move(evtv);
}

void DDGev::DhCP(uint64_t ncpp, uint64_t ncpn, const AbsDDPars& pars,
                 const string& fname) {
    auto evtv = DhCP(ncpp, pars, dtypes::DhCPp);
    auto evtcpn = DhCP(ncpn, pars, dtypes::DhCPn);
    evtv->insert(evtv->end(),
                 make_move_iterator(evtcpn->begin()),
                 make_move_iterator(evtcpn->end()));
    cout << evtv->size() << " / " << ncpp + ncpn
         << " events generated" << endl;
    // Random shuffle positive and negative tags
    std::random_shuffle(evtv->begin(), evtv->end());
    // Write to a text file
    ofstream ofile(fname, ofstream::out);
    for (const auto& evt : *evtv)
        ofile << *evt << endl;
    ofile.close();
}

void DDGev::DhDalitz(uint64_t nevt, const AbsDDPars& pars,
                     const string& fname) {
    auto pdf = Cfg::pdf(m_expCfg);  // Toy PDF
    // Decay rates //
    double rate_norm = 0;
    DMap rates;
    const double alpha = pdf->alpha();
    for (auto dbin = -8; dbin <= 8; dbin++) if (dbin) {
        auto ud = pars.rawud(0, dbin, dtypes::Dh);
        auto rate = ud.first + alpha * ud.second;
        rate_norm += rate;
        rates[dbin] = rate;
    }
    cout << "DhDalitz: " << endl
         << "  alpha: " << alpha << endl
         << "  rnorm: " << rate_norm << endl;
    libTatami::ToyPdfGen gen(*pdf);  // generator
    vector<unique_ptr<BEvt>> evtv;  // vector of events
    for (auto dbin = -8; dbin <= 8; dbin++) if (dbin) {
        cout << "Bin " << dbin << endl;
        auto cs = pars.coefs(0, dbin, dtypes::Dh);
        if (!isgood_cs(cs)) {
            cerr << "DDGev::DhDalitz Dbin: " << dbin << endl;
            throw new std::runtime_error("Bad C or S");
        }
        if (dump) cout << "c " << cs.first << ", s " << cs.second << endl;
        pdf->SetCS(cs);
        auto half_nevt = static_cast<uint64_t>
                (0.5 * nevt * rates[dbin] / rate_norm);
        pdf->SetTag(1);
        if (dump) cout << "DDGev::DhDalitz: time generation" << endl;
        auto times = gen.Generate(half_nevt, true);
        if (dump) cout << "DDGev::DhDalitz: time generated " << times.size() << endl;
        for (auto t : times)
            evtv.emplace_back(make_unique<BEvt>(t, 1, dbin, dtypes::Dh));
        pdf->SetTag(-1);
        times = gen.Generate(half_nevt, true);
        for (auto t : times)
            evtv.emplace_back(make_unique<BEvt>(t, -1, dbin, dtypes::Dh));
    }
    cout << evtv.size() << " / " << nevt << " events generated" << endl;
    // Random shuffle positive and negative tags
    std::random_shuffle(evtv.begin(), evtv.end());
    // Write to a text file
    ofstream ofile(fname, ofstream::out);
    for (const auto& evt : evtv)
        ofile << *evt << endl;
    ofile.close();
}

void DDGev::DoubleDalitz(uint64_t nevt, const AbsDDPars& pars,
                          const string& fname) {
    auto pdf = Cfg::pdf(m_expCfg);  // Toy PDF
    // Decay rates //
    double rate_norm = 0;
    DDMap rates;
    const double alpha = pdf->alpha();
    for (auto bbin = -8; bbin <= 8; bbin++) if (bbin) {
        for (auto dbin = -8; dbin <= 8; dbin++) if (dbin) {
            auto ud = pars.rawud(bbin, dbin, dtypes::KsPIPI);
            auto rate = ud.first + alpha * ud.second;
            rate_norm += rate;
            rates[bbin][dbin] = rate;
        }
    }
    cout << "DoubleDalitzWithWF: " << endl
         << "  alpha: " << alpha << endl
         << "  rnorm: " << rate_norm << endl;
    libTatami::ToyPdfGen gen(*pdf);  // generator
    vector<unique_ptr<DBEvt>> evtv;  // vector of events
    for (auto bbin = -8; bbin <= 8; bbin++) if (bbin) {
        for (auto dbin = -8; dbin <= 8; dbin++) if (dbin) {
            auto cs = pars.coefs(bbin, dbin, dtypes::KsPIPI);
            if (!isgood_cs(cs)) {
                cerr << "DDGev::DoubleDalitz Bbin: " << bbin
                     << ", Dbin: " << dbin << endl;
                throw new std::runtime_error("Bad C or S");
            }
            pdf->SetCS(cs);
            auto half_nevt = static_cast<uint64_t>
                    (0.5 * nevt * rates[bbin][dbin] / rate_norm);
            pdf->SetTag(1);
            auto times = gen.Generate(half_nevt, true);
            for (auto t : times)
                evtv.emplace_back(make_unique<DBEvt>(t, 1, dbin, bbin));
            pdf->SetTag(-1);
            times = gen.Generate(half_nevt, true);
            for (auto t : times)
                evtv.emplace_back(make_unique<DBEvt>(t, -1, dbin, bbin));
        }
    }
    cout << evtv.size() << " / " << nevt << " events generated" << endl;
    // Random shuffle positive and negative tags
    std::random_shuffle(evtv.begin(), evtv.end());
    // Write to a text file
    ofstream ofile(fname, ofstream::out);
    for (const auto& evt : evtv)
        ofile << *evt << endl;
    ofile.close();
}
