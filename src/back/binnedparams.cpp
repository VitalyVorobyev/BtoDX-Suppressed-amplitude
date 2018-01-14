#include "binnedparams.h"

#include <utility>      // std::forward
#include <iostream>
#include <cmath>
#include <iomanip>

using std::sqrt;
using std::sin;
using std::cos;

using std::vector;
using std::string;

using std::cout;
using std::cin;
using std::endl;

using std::make_pair;
using std::pair;

using vectd = vector<double>;

BinnedParams::BinnedParams(uint16_t n) :
    nbins(n), ccoefpcl(2*nbins), scoefpcl(2*nbins), updated(false) {
    vector<string> names{"Kp", "Kn", "C", "S", "Kt", "C1", "S1", "C2", "S2"};
    for (string& name : names)
        params.emplace(make_pair(name, vectd(nbins)));
}

auto BinnedParams::prc(std::string&& label) const {
    return precl.find(label)->second;
}

auto BinnedParams::get(string&& label, int16_t bin) const {
    return params.find(label)->second[bin];
}

auto BinnedParams::ccoef(double Kp, double Kn, double C, double S,
                         double C1, double S1) const {
    double cmb1 = S  * prc("cosgam") - C  * prc("singam");
    double cmb2 = C1 * prc("cosgam") + S1 * prc("singam");
    double inbr = prc("sindb") * cmb1 - prc("cosdb") * cmb2;
    return (Kp - Kn + 4. * rb * sqrt(Kp * Kn) * inbr) / (Kp + Kn);
}

auto BinnedParams::scoef(double Kp, double Kn, double Kt,
                         double C, double S, double C2, double S2) const {
    double sqkk = sqrt(Kp * Kn);
    double cmb1 = S * prc("cos2beta") + C * prc("sin2beta");
    double cmb2 = Kp * prc("sin_g+2b-db") + Kn * prc("sin_g+2b+db");
    double cmb3 = Kt * prc("sin_g+2b") + S2 * prc("cos_2b-g")
                                       + C2 * prc("sin_2b-g");
    return 2. * (-sqkk*cmb1 - rb*cmb2 + 2.*rb*sqkk*prc("cosdb")*cmb3)
            / (Kp + Kn);
}

void BinnedParams::update() {
    // Trigonometry caching
    precl["singam"] = sin(gamma);
    precl["cosgam"] = cos(gamma);
    precl["cosdb"] = cos(delb);
    precl["sindb"] = sin(delb);
    precl["cos2beta"] = cos(2.*beta);
    precl["sin2beta"] = sin(2.*beta);
    precl["sin_g+2b-db"] = sin(gamma + 2.*beta - delb);
    precl["sin_g+2b+db"] = sin(gamma + 2.*beta + delb);
    precl["sin_g+2b"] = sin(2.*beta + gamma);
    precl["sin_2b-g"] = sin(2.*beta - gamma);
    precl["cos_2b-g"] = cos(2.*beta - gamma);

    for (auto bin = 0; bin < nbins; bin++) {
        // positive bin indexes
        ccoefpcl[bin+1] = ccoef(get("Kp", bin), get("Kn", bin),
                get("C", bin), get("S", bin), get("C1", bin), get("S1", bin));
        scoefpcl[bin+1] = scoef(get("Kp", bin), get("Kn", bin), get("Kt", bin),
                get("C", bin), get("S", bin), get("C2", bin), get("S2", bin));
        // negative bin indexes
        ccoefpcl[-bin-1] = ccoef(get("Kn", bin), get("Kp", bin),
                get("C", bin), -get("S", bin), -get("C1", bin), get("S1", bin));
        scoefpcl[-bin-1] = scoef(get("Kn", bin), get("Kp", bin), get("Kt", bin),
                get("C", bin), -get("S", bin), get("C2", bin), -get("S2", bin));
    }
    updated = true;
}

void BinnedParams::setBeta(double b) {
    beta = b;
    updated = false;
}

void BinnedParams::setGamma(double b) {
    gamma = b;
    updated = false;
}

void BinnedParams::setRb(double b) {
    rb = b;
    updated = false;
}

void BinnedParams::setDeltaB(double b) {
    delb = b;
    updated = false;
}

void BinnedParams::setParam(string label, uint16_t bin, double val) {
    if (params.count(label) && bin && (bin < nbins)) {
        params.find(label)->second[bin-1] = val;
        updated = false;
    } else {
        cout << "Bad setParam call: " << label
             << " " << bin << " " << val << endl;
    }
}

void BinnedParams::setParam(string label, vectd parv) {
    params[label] = parv;
    updated = false;
}

pair<double, double> BinnedParams::coefs(int16_t bin) {
    // bin is in {-nbins, -nbins+1, ..., -1, 1, 2, ..., nbins}
    if (!bin || (abs(bin) > nbins)) {
        cout << "BinnedParams::coefs: wrong index " << bin << "/" << nbins << endl;
        return make_pair(0., 0.);
    }
    if (!updated) update();
    return make_pair(ccoefpcl[bin], scoefpcl[bin]);
}

pair<double, double> BinnedParams::EvtFrac(int16_t bin) const {
    if (!bin || (abs(bin) > nbins)) {
        cout << "BinnedParams::coefs: wrong index " << bin << "/" << nbins << endl;
        return make_pair(0., 0.);
    }
    if (bin < 0)
        return make_pair(get("Kn", -bin - 1), get("Kp", -bin - 1));
    return make_pair(get("Kp", bin - 1), get("Kn", bin - 1));
}

void BinnedParams::print() {
    if (!updated) update();
    cout << "### Binned coeffs list ###" << endl;
    for (auto bin = 1; bin <= 8; bin++) {
        cout << "  " << bin << " ->"
             << " c+: " << std::setw(5) << ccoefpcl[bin]
             << " s+: " << std::setw(5) << scoefpcl[bin]
             << " c-: " << std::setw(5) << ccoefpcl[-bin]
             << " s-: " << std::setw(5) << scoefpcl[-bin]
             << endl;
    }
}
