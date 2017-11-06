#include "cpvgen.h"

#include <algorithm>
#include <iostream>
#include <fstream>
#include <map>

#include "mylibs/libTatami/toypdfgen.h"

#include "cfg.h"

constexpr bool dump = false;

const std::string CPVGen::default_str("");

using std::ofstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;

using std::move;
using std::unique_ptr;
using std::make_unique;
using std::make_move_iterator;

auto CPVGen::Gen(libTatami::ToyPdf& pdf, pair<uint64_t, uint64_t>&& nev,
         dtypes type, int16_t bin) {
    vector<unique_ptr<BEvt>> evtv;

    libTatami::ToyPdfGen gen(pdf);
    // Generate dt distributions
    pdf.SetTag(1);
    auto timev = gen.Generate(nev.first, true);
    for (auto time : timev)
        evtv.emplace_back(make_unique<BEvt>(time, +1, bin, type));
    pdf.SetTag(-1);
    timev = gen.Generate(nev.second, true);
    for (auto time : timev)
        evtv.emplace_back(make_unique<BEvt>(time, -1, bin, type));
    return move(evtv);
}

auto GenTime(const libTatami::ToyPdf& pdf, uint64_t nevt) {
    libTatami::ToyPdfGen gen(pdf);
    return gen.Generate(nevt, true);
}

auto CPVGen::Gen(dtypes type) {
    auto pdf = Cfg::pdf();
    auto lamf = Cfg::lamf(type);
    auto nev = Cfg::nevts(type);
    pdf.SetS(lamf.scoef());
    pdf.SetC(lamf.ccoef());
    return Gen(pdf, nev, type);
}

void CPVGen::ShuffleAndWrite(vector<unique_ptr<BEvt>>& evts,
                             const string& fname) {
    // Random shuffle positive and negative tags
    std::random_shuffle(evts.begin(), evts.end());

    // Write to text
    ofstream ofile(fname, ofstream::out);
    for (const auto& evt : evts)
        ofile << evt->t() << " " << evt->Tag() << " " << evt->Bin() << " "
              << static_cast<int>(evt->Type()) << endl;
    ofile.close();
}

void CPVGen::GenAndWrite(libTatami::ToyPdf& pdf,
                         pair<uint64_t, uint64_t> &&nev,
                         const string& fname, dtypes type) {
    auto evtv = Gen(pdf, move(nev), type);
    ShuffleAndWrite(evtv, fname);
}

void CPVGen::GenAndWrite(dtypes type, const string& pref) {
    auto evtv = Gen(type);
    auto fname = Cfg::dfile(type) + pref;
    ShuffleAndWrite(evtv, fname);
}

void CPVGen::GenAndWrite(dtypes type, pair<uint64_t, uint64_t>&& nev,
                         const string &fname) {
    if (type == dtypes::KsPIPI) {
        GenKspipi(move(nev), fname);
        return;
    }
    auto pdf = Cfg::pdf();
    auto lamf = Cfg::lamf(type);
    pdf.SetS(lamf.scoef());
    pdf.SetC(lamf.ccoef());
    auto filename = (fname == default_str) ? Cfg::dfile(type, nev.first) :
                                             fname;
    cout << filename << endl << Cfg::dfile(type, nev.first) << endl;
    GenAndWrite(pdf, move(nev), filename, type);
}

void CPVGen::GenKspipi(pair<uint64_t, uint64_t>&& nev, const string& fname) {
    const auto bpars = Cfg::bpars();
    bpars->print();
    vector<unique_ptr<BEvt>> evtv;
    auto pdf = Cfg::pdf();
    for (auto bin = 1; bin <= 8; bin++) {
        const auto fracs = bpars->EvtFrac(bin);

        // positive bin
        const auto coefsp = bpars->coefs(bin);
        pdf.SetC(coefsp.first);
        pdf.SetS(coefsp.second);
        auto evtnp = make_pair(nev.first*fracs.first, nev.second*fracs.first);
        auto evtvec = Gen(pdf, evtnp, dtypes::KsPIPI, bin);
        evtv.insert(evtv.end(), make_move_iterator(evtvec.begin()),
                                make_move_iterator(evtvec.end()));
        // negative bin
        const auto coefsn = bpars->coefs(-bin);
        pdf.SetC(coefsn.first);
        pdf.SetS(coefsn.second);
        auto evtnn = make_pair(nev.first*fracs.second,
                               nev.second*fracs.second);
        evtvec = Gen(pdf, evtnn, dtypes::KsPIPI, -bin);
        evtv.insert(evtv.end(), make_move_iterator(evtvec.begin()),
                                make_move_iterator(evtvec.end()));
    }
    cout << evtv.size() << " events genetated" << endl;
    if (fname == default_str) {
        cout << "Write events in "
             << Cfg::dfile(dtypes::KsPIPI, nev.first) << endl;
        ShuffleAndWrite(evtv, Cfg::dfile(dtypes::KsPIPI, nev.first));
    } else {
        cout << "Write events in " << fname << endl;
        ShuffleAndWrite(evtv, fname);
    }
}

void CPVGen::CompleteData(const string& pref) {
    GenAndWrite(dtypes::CPp, pref);
    GenAndWrite(dtypes::CPn, pref);
    GenAndWrite(dtypes::KPI, pref);
    GenAndWrite(dtypes::PIK, pref);
}
