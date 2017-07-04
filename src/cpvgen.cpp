#include "cpvgen.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "mylibs/libTatami/toypdfgen.h"

#include "cfg.h"

const std::string CPVGen::default_str("");

using std::ofstream;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::make_pair;
using std::move;

void CPVGen::Gen(libTatami::ToyPdf& pdf, pair<uint64_t, uint64_t>&& nev,
         vector<BEvt>& evtv, dtypes type, int16_t bin) {
    libTatami::ToyPdfGen gen(&pdf);
//    evtv.reserve(nev.first + nev.second);
    std::vector<double> timev;
    // Generate dt distributions
    pdf.SetTag(1);
    gen.Generate(nev.first, &timev, true);
    for (auto time : timev) evtv.emplace_back(BEvt(time, +1, bin, type));
    pdf.SetTag(-1);
    timev.clear();
    gen.Generate(nev.second, &timev, true);
    for (auto time : timev) evtv.emplace_back(BEvt(time, -1, bin, type));
}

void CPVGen::Gen(vector<BEvt>& evtv, dtypes type) {
    auto pdf = Cfg::pdf();
    auto lamf = Cfg::lamf(type);
    auto nev = Cfg::nevts(type);
    pdf.SetS(lamf.scoef());
    pdf.SetC(lamf.ccoef());
    Gen(pdf, nev, evtv, type);
}

void CPVGen::ShuffleAndWrite(vector<BEvt>& evts, const string& fname) {
    // Random shuffle positive and negative tags
    std::random_shuffle(evts.begin(), evts.end());

    // Write to text
    ofstream ofile(fname, ofstream::out);
    for (auto& evt : evts)
        ofile << evt.time << " " << evt.tag << " " << evt.bin << " "
              << static_cast<int>(evt.type) << endl;
    ofile.close();
}

void CPVGen::GenAndWrite(libTatami::ToyPdf& pdf,
                         pair<uint64_t, uint64_t> &&nev,
                         const string& fname, dtypes type) {
    std::vector<BEvt> evtv;
    Gen(pdf, move(nev), evtv, type);
    ShuffleAndWrite(evtv, fname);
}

void CPVGen::GenAndWrite(dtypes type, const string& pref) {
    std::vector<BEvt> evtv;
    Gen(evtv, type);
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
    auto filename = (fname == default_str) ? Cfg::dfile(type, nev.first) : fname;
    cout << filename << endl << Cfg::dfile(type, nev.first) << endl;
    GenAndWrite(pdf, move(nev), filename, type);
}

void CPVGen::GenKspipi(pair<uint64_t, uint64_t>&& nev, const string& fname) {
    const auto bpars = Cfg::bpars();
    bpars->print();
    vector<BEvt> evtv;
    auto pdf = Cfg::pdf();
    for (auto bin = 1; bin <= 8; bin++) {
        const auto fracs = bpars->EvtFrac(bin);

        // positive bin
        const auto coefsp = bpars->coefs(bin);
        pdf.SetC(coefsp.first);
        pdf.SetS(coefsp.second);
        auto evtnp = make_pair(nev.first*fracs.first, nev.second*fracs.first);
        Gen(pdf, evtnp, evtv, dtypes::KsPIPI, bin);

        // negative bin
        const auto coefsn = bpars->coefs(-bin);
        pdf.SetC(coefsn.first);
        pdf.SetS(coefsn.second);
        auto evtnn = make_pair(nev.first*fracs.second, nev.second*fracs.second);
        Gen(pdf, evtnn, evtv, dtypes::KsPIPI, -bin);
    }
    cout << evtv.size() << " events genetated" << endl;
    if (fname == default_str) {
        cout << "Write events in " << Cfg::dfile(dtypes::KsPIPI, nev.first) << endl;
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
