#include "absddpars.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>

using std::string;
using std::vector;
using std::ifstream;
using std::cout;
using std::cerr;
using std::endl;
using std::setfill;
using std::setw;
using std::setprecision;

AbsDDPars::AbsDDPars(const string& dcfg, const string& bcfg) {
    readDConfig(dcfg, true);
    readBConfig(bcfg, true);
}

void AbsDDPars::set_c(vector<double>& x) {
    m_int["Crf"] = std::move(x);
}

void AbsDDPars::set_s(vector<double>& x) {
    m_int["Srf"] = std::move(x);
}

int16_t AbsDDPars::readDConfig(const string& fname, bool verb) {
    if (verb) cout << "### AbsDDPars::readDConfig ###" << endl
                   << fname << endl;
    ifstream cfg(fname, ifstream::in);
    if (!cfg.good()) {
        cerr << "Can't open file " << fname << endl;
        return -1;
    }
    vector<string> labels{"K+", "K-", "C", "S"};
    for (auto& lbl : labels)
        m_int.emplace(lbl, vector<double>(ndbins()));

    for (auto idx = 0u; idx < ndbins(); idx++) {
        static string line;
        static uint32_t bin;
        static double C, S, Kp, Kn;

        getline(cfg, line);
        auto flag = sscanf(line.c_str(),
                    "%u: C = %lf, S = %lf, K+ = %lf, K- = %lf",
                    &bin, &C, &S, &Kp, &Kn);
        if (flag != 5) return -2;
        m_int["C"][idx] = C;
        m_int["S"][idx] = S;
        m_int["K+"][idx] = Kp;
        m_int["K-"][idx] = Kn;
    }
    return 0;
}

int16_t AbsDDPars::readBConfig(const string& fname, bool verb) {
    if (verb) cout << "### AbsDDPars::readBConfig ###" << endl
                   << fname << endl;
    ifstream cfg(fname, ifstream::in);
    if (!cfg.good()) {
        cerr << "Can't open file " << fname << endl;
        return -1;
    }
    vector<string> labels{"K+rf", "K-rf", "K+wf", "K-wf", "Crf", "Srf", "Cwf", "Swf",
                       "Ctp", "Stp", "Ctn", "Stn", "Cpp", "Spp", "Cpn", "Spn"};
    for (auto& lbl : labels)
        m_int.emplace(lbl, vector<double>(nbbins()));

    for (auto idx = 0u; idx < nbbins(); idx++) {
        static string line;
        static uint32_t bin;
        static double p1, p2, p3, p4;

        getline(cfg, line);
        auto flag = sscanf(line.c_str(),
                    "%u: K+rf = %lf, K-rf = %lf, K+wf = %lf, K-wf = %lf",
                    &bin, &p1, &p2, &p3, &p4);
        if (flag != 5) return -2;
        m_int["K+rf"][idx] = p1;
        m_int["K-rf"][idx] = p2;
        m_int["K+wf"][idx] = p3;
        m_int["K-wf"][idx] = p4;

        getline(cfg, line);
        flag = sscanf(line.c_str(),
                      "  Crf = %lf, Srf = %lf, Cwf = %lf, Swf = %lf",
                      &p1, &p2, &p3, &p4);
        if (flag != 4) return -3;
        m_int["Crf"][idx] = p1;
        m_int["Srf"][idx] = p2;
        m_int["Cwf"][idx] = p3;
        m_int["Swf"][idx] = p4;

        getline(cfg, line);
        flag = sscanf(line.c_str(),
                      "  Ctp = %lf, Stp = %lf, Ctn = %lf, Stn = %lf",
                      &p1, &p2, &p3, &p4);
        if (flag != 4) return -4;
        m_int["Ctp"][idx] = p1;
        m_int["Stp"][idx] = p2;
        m_int["Ctn"][idx] = p3;
        m_int["Stn"][idx] = p4;

        getline(cfg, line);
        flag = sscanf(line.c_str(),
                      "  Cpp = %lf, Spp = %lf, Cpn = %lf, Spn = %lf",
                      &p1, &p2, &p3, &p4);
        if (flag != 4) return -5;
        m_int["Cpp"][idx] = p1;
        m_int["Spp"][idx] = p2;
        m_int["Cpn"][idx] = p3;
        m_int["Spn"][idx] = p4;
    }
    return 0;
}

inline void formOut(const string& s, int pr, int wd) {
    cout << setw(wd) << setprecision(pr) << setfill(' ') << s;
}

void AbsDDPars::print() const {
    cout << "### AbsDDPars::DConfig ###" << endl;
    for (auto idx = 0u; idx < ndbins(); idx++) {
        cout << "  " << idx + 1 << ": "
             << "C = " << setw(5) << setfill(' ') << m_int.at("C")[idx] << ", "
             << "S = " << setw(5) << setfill(' ') << m_int.at("S")[idx] << ", "
             << "K+ = " << setw(5) << setfill(' ') << m_int.at("K+")[idx] << ", "
             << "K- = " << setw(5) << setfill(' ') << m_int.at("K-")[idx] << endl;
    }
    cout << endl << "### AbsDDPars::BConfig ###" << endl;
    auto pr = 6;
    auto wd = 9;
    for (auto idx = 0u; idx < nbbins(); idx++) {
        cout << "    " << idx + 1 << ": "
             << "K+rf = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("K+rf")[idx] << ", "
             << "K-rf = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("K-rf")[idx] << ", "
             << "K+wf = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("K+wf")[idx] << ", "
             << "K-wf = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("K-wf")[idx] << endl;
        cout << "       Crf  = " << std::setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Crf")[idx] << ", "
             << "Srf  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Srf")[idx] << ", "
             << "Cwf  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Cwf")[idx] << ", "
             << "Swf  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Swf")[idx] << endl;
        cout << "       Ctp  = " << std::setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Ctp")[idx] << ", "
             << "Stp  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Stp")[idx] << ", "
             << "Ctn  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Ctn")[idx] << ", "
             << "Stn  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Stn")[idx] << endl;
        cout << "       Cpp  = " << std::setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Cpp")[idx] << ", "
             << "Spp  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Spp")[idx] << ", "
             << "Cpn  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Cpn")[idx] << ", "
             << "Spn  = " << setw(wd) << setprecision(pr) << setfill(' ') << m_int.at("Spn")[idx] << endl;
    }
}

const std::vector<double>& AbsDDPars::get_c() const {
    if (m_int.find("Crf") == m_int.end())
        throw new std::runtime_error("get_c: Crf is not in m_int");
    return m_int.at("Crf");
}
const std::vector<double>& AbsDDPars::get_s() const {
    if (m_int.find("Srf") == m_int.end())
        throw new std::runtime_error("get_s: Srf is not in m_int");
    return m_int.at("Srf");
}
