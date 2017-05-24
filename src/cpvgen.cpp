#include "cpvgen.h"

#include <algorithm>
#include <iostream>
#include <fstream>

#include "mylibs/libTatami/toypdfgen.h"

void CPVGen::Gen(libTatami::ToyPdf& pdf, std::pair<int, int> nev,
         std::vector<Evt>& evtv, dtypes type) {
    libTatami::ToyPdfGen gen(&pdf);
    evtv.reserve(nev.first + nev.second);
    std::vector<double> timev;
    // Generate dt distributions
    pdf.SetTag(1);
    gen.Generate(nev.first, &timev, true);
    for (auto time : timev) evtv.push_back(Evt(time, +1, type));
    pdf.SetTag(-1);
    timev.clear();
    gen.Generate(nev.second, &timev, true);
    for (auto time : timev) evtv.push_back(Evt(time, -1, type));
}

void CPVGen::Gen(std::vector<Evt>& evtv, dtypes type) {
    auto pdf = Cfg::pdf();
    auto lamf = Cfg::lamf(type);
    auto nev = Cfg::nevts(type);
    pdf.SetS(lamf.scoef());
    pdf.SetC(lamf.ccoef());
    Gen(pdf, nev, evtv, type);
}

void CPVGen::ShuffleAndWrite(std::vector<Evt>& evts, std::string& fname) {
    // Random shuffle positive and negative tags
    std::random_shuffle(evts.begin(), evts.end());

    // Write to text
    std::ofstream ofile(fname, std::ofstream::out);
    for (auto evt : evts)
        ofile << evt.time << " " << evt.tag << " " << evt.type << std::endl;
    ofile.close();
}

void CPVGen::GenAndWrite(libTatami::ToyPdf& pdf, std::pair<int, int> nev,
                 std::string& fname, dtypes type) {
    std::vector<Evt> evtv;
    Gen(pdf, nev, evtv, type);
    ShuffleAndWrite(evtv, fname);
}

void CPVGen::GenAndWrite(dtypes type) {
    std::vector<Evt> evtv;
    Gen(evtv, type);
    auto fname = Cfg::dfile(type);
    ShuffleAndWrite(evtv, fname);
}

void CPVGen::GenAndWrite(dtypes type, std::pair<int, int> nev) {
    std::vector<Evt> evtv;
    auto pdf = Cfg::pdf();
    auto fname = Cfg::dfile(type);
    auto lamf = Cfg::lamf(type);
    pdf.SetS(lamf.scoef());
    pdf.SetC(lamf.ccoef());
    GenAndWrite(pdf, nev, fname, type);
}

void CPVGen::CompleteData() {
    GenAndWrite(dtypes::CPp);
    GenAndWrite(dtypes::CPn);
    GenAndWrite(dtypes::KPI);
    GenAndWrite(dtypes::PIK);
}
