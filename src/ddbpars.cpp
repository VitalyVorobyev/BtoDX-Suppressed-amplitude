#include "ddbpars.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using str = std::string;
using vectd = std::vector<double>;

using std::ifstream;

using std::cout;
using std::endl;
using std::cerr;

using std::sin;
using std::cos;
using std::sqrt;
using std::getline;
using std::sscanf;
using std::make_pair;
using std::vector;

constexpr bool dump = false;
constexpr double epsilon = 0.00001;

DDBPars::DDBPars(double rb, double beta, double gamma, double delb,
                 const str& dcfg, const str& bcfg) : AbsDDPars(dcfg, bcfg) {
    set(rb, beta, gamma);
    set_delb(delb);
    print();
}

void DDBPars::set(double rb, double _beta, double gamma) {
    m_rb = rb;
    AbsDDPars::set_beta(_beta);
    m_gamma = gamma;
    cout << "DDBPars: rb = " << m_rb
         << ", beta = " << beta()
         << ", gamma = " << m_gamma << endl;
    m_cache["rb2"] = m_rb * m_rb;
    m_cache["cos(2beta)"] = cos(2. * beta());
    m_cache["sin(2beta)"] = sin(2. * beta());
    m_cache["cos(gamma)"] = cos(m_gamma);
    m_cache["sin(gamma)"] = sin(m_gamma);
    cache_upd();
}

void DDBPars::set_rb(double x) {
    m_rb = x;
    m_cache["rb2"] = m_rb * m_rb;
}

void DDBPars::set_delb(double x) {
    cout << "DDBPars delb = " << x << endl;
    m_delb = x;
    m_cache["sin(delb)"] = sin(x);
    m_cache["cos(delb)"] = cos(x);
}

void DDBPars::cache_upd() {
    m_cache["cos(2beta + gamma)"] = cos(2. * beta() + m_gamma);
    m_cache["sin(2beta + gamma)"] = sin(2. * beta() + m_gamma);
    m_cache["cos(2beta + 2gamma)"] = cos(2. * (beta() + m_gamma));
    m_cache["sin(2beta + 2gamma)"] = sin(2. * (beta() + m_gamma));
}

void DDBPars::set_beta(double x) {
    AbsDDPars::set_beta(x);
    m_cache["cos(2beta)"] = cos(2. * beta());
    m_cache["sin(2beta)"] = sin(2. * beta());
    cache_upd();
}

void DDBPars::set_gamma(double x) {
    m_gamma = x;
    m_cache["cos(gamma)"] = cos(m_gamma);
    m_cache["sin(gamma)"] = sin(m_gamma);
    cache_upd();
}

void DDBPars::set_pars(int16_t bbin, int16_t dbin) const {
    if (!bbin && !dbin) return;
    if (dump) cout << "Get pars" << endl;
    // D meson binned Dalitz params
    if (dbin > 0) {
        dbin--;  // bin to index
        m_pars.Kp = m_int.at("K+")[dbin];
        m_pars.Kn = m_int.at("K-")[dbin];
        m_pars.C = m_int.at("C")[dbin];
        m_pars.S = m_int.at("S")[dbin];
    } else if (dbin < 0) {
        dbin = -dbin - 1;  // bin to index
        m_pars.Kp = m_int.at("K-")[dbin];
        m_pars.Kn = m_int.at("K+")[dbin];
        m_pars.C = m_int.at("C")[dbin];
        m_pars.S = -m_int.at("S")[dbin];
    }
    if (bbin > 0) {
        bbin--;
        m_pars.Kprf = m_int.at("K+rf")[bbin];
        m_pars.Knrf = m_int.at("K-rf")[bbin];
        m_pars.Crf = m_int.at("Crf")[bbin];
        m_pars.Srf = m_int.at("Srf")[bbin];
    } else if (bbin < 0) {
        bbin = -bbin - 1;
        m_pars.Kprf = m_int.at("K-rf")[bbin];
        m_pars.Knrf = m_int.at("K+rf")[bbin];
        m_pars.Crf = m_int.at("Crf")[bbin];
        m_pars.Srf = -m_int.at("Srf")[bbin];
    }
    if (m_rb < epsilon) return;
    // B meson binned Dalitz params
    if (bbin > 0) {
        m_pars.Kpwf = m_int.at("K+wf")[bbin];
        m_pars.Knwf = m_int.at("K-wf")[bbin];
        m_pars.Cwf = m_int.at("Cwf")[bbin];
        m_pars.Swf = m_int.at("Swf")[bbin];

        m_pars.Ctp = m_int.at("Ctp")[bbin];
        m_pars.Stp = m_int.at("Stp")[bbin];
        m_pars.Ctn = m_int.at("Ctn")[bbin];
        m_pars.Stn = m_int.at("Stn")[bbin];

        m_pars.Cpp = m_int.at("Cpp")[bbin];
        m_pars.Spp = m_int.at("Spp")[bbin];
        m_pars.Cpn = m_int.at("Cpn")[bbin];
        m_pars.Spn = m_int.at("Spn")[bbin];
    } else if (bbin < 0) {
        m_pars.Kpwf = m_int.at("K-wf")[bbin];
        m_pars.Knwf = m_int.at("K+wf")[bbin];
        m_pars.Cwf = m_int.at("Cwf")[bbin];
        m_pars.Swf = -m_int.at("Swf")[bbin];

        m_pars.Ctp = m_int.at("Ctn")[bbin];
        m_pars.Stp = m_int.at("Stn")[bbin];
        m_pars.Ctn = m_int.at("Ctp")[bbin];
        m_pars.Stn = m_int.at("Stp")[bbin];

        m_pars.Cpp = m_int.at("Cpn")[bbin];
        m_pars.Spp = m_int.at("Spn")[bbin];
        m_pars.Cpn = m_int.at("Cpp")[bbin];
        m_pars.Spn = m_int.at("Spp")[bbin];
    }
}

DDBPars::ddpair DDBPars::calc_ud() const {
    if (dump) cout << "calc_ud" << endl;
    // Calculate U //
    auto KpKn = m_pars.Kp * m_pars.Kn;
    auto p11 = m_pars.Kn * m_pars.Kprf;
    auto p12 = m_pars.Kp * m_pars.Knrf;
    auto p1 = 0.5 * (p11 + p12);

    if (m_rb < epsilon)
        return make_pair(p1, 0.5 * (p11 - p12));

    auto p20 = m_rb * sqrt(KpKn * m_pars.Kprf * m_pars.Kpwf);
    auto p21 = (m_pars.Ctp * m_pars.C +
                m_pars.Stp * m_pars.S) * m_cache.at("cos(gamma)");
    auto p22 = (m_pars.Ctp * m_pars.S +
                m_pars.Stp * m_pars.C) * m_cache.at("sin(gamma)");
    auto p2 = p20 * (p21 - p22);

    auto p30 = m_rb * sqrt(KpKn * m_pars.Knrf * m_pars.Knwf);
    auto p31 = (m_pars.Ctn * m_pars.C -
                m_pars.Stn * m_pars.S) * m_cache.at("cos(gamma)");
    auto p32 = (m_pars.Ctn * m_pars.S +
                m_pars.Stn * m_pars.C) * m_cache.at("sin(gamma)");
    auto p3 = p30 * (p31 - p32);

    auto p41 = m_pars.Kp * m_pars.Kpwf;
    auto p42 = m_pars.Kn * m_pars.Knwf;
    auto p4 = 0.5 * m_cache.at("rb2") * (p41 + p42);

    auto u = p1 + p2 + p3 + p4;

    // Calcuate D //
    p1 = 0.5 * (p11 - p12);
    p4 = 0.5 * m_cache.at("rb2") * (p41 - p42);
    auto d = p1 + p2 - p3 + p4;
    return make_pair(u, d);
}

double DDBPars::calc_f() const {
    if (dump) cout << "calc_f" << endl;
    // Calculate F //
    auto KpKn = m_pars.Kp * m_pars.Kn;
    auto p10 = sqrt(KpKn * m_pars.Kprf * m_pars.Knrf);
    auto p11 = (m_pars.Srf * m_pars.C -
                m_pars.Crf * m_pars.S) * m_cache.at("cos(2beta)");
    auto p12 = (m_pars.Crf * m_pars.C +
                m_pars.Srf * m_pars.S) * m_cache.at("sin(2beta)");
    auto p1 = -p10 * (p11 - p12);
    if (m_rb < epsilon) return p1;

    auto p20 = m_rb * m_pars.Kn * sqrt(m_pars.Kprf * m_pars.Knwf);
    auto p21 = m_pars.Cpp * m_cache.at("sin(2beta + gamma)");
    auto p22 = m_pars.Spp * m_cache.at("cos(2beta + gamma)");
    auto p2 = p20 * (p21 - p22);

    auto p30 = m_rb * m_pars.Kp * sqrt(m_pars.Knrf * m_pars.Kpwf);
    auto p31 = m_pars.Cpn * m_cache.at("sin(2beta + gamma)");
    auto p32 = m_pars.Spn * m_cache.at("cos(2beta + gamma)");
    auto p3 = p30 * (p31 + p32);

    auto p40 = 0.5 * m_cache.at("rb2") *
            sqrt(KpKn * m_pars.Kpwf * m_pars.Knwf);
    auto p41 = (m_pars.Swf * m_pars.C +
                m_pars.Cwf * m_pars.S) * m_cache.at("cos(2beta + 2gamma)");
    auto p42 = (m_pars.Cwf * m_pars.C -
                m_pars.Swf * m_pars.S) * m_cache.at("sin(2beta + 2gamma)");
    auto p4 = -p40 * (p41 - p42);
    return p1 + p2 + p3 + p4;
}

DDBPars::ddpair DDBPars::calc_ud_cp(int xiD) const {
    if (dump) cout << "calc_ud_cp" << endl;
    // Calculate U and D //
    auto Xplus = 0.5 * (m_pars.Kprf + m_cache.at("rb2") * m_pars.Kpwf);
    auto Xmins = 0.5 * (m_pars.Knrf + m_cache.at("rb2") * m_pars.Knwf);
    if (m_rb >= epsilon) {
        Xplus += xiD * m_rb * sqrt(m_pars.Kprf * m_pars.Kpwf) *
                (m_pars.Ctp * m_cache.at("cos(gamma)") -
                 m_pars.Stp * m_cache.at("sin(gamma)"));
        Xmins += xiD * m_rb * sqrt(m_pars.Knrf * m_pars.Knwf) *
                (m_pars.Ctn * m_cache.at("cos(gamma)") +
                 m_pars.Stn * m_cache.at("sin(gamma)"));
    }
    return make_pair(Xplus + Xmins, Xplus - Xmins);
}

double DDBPars::calc_f_cp(int xiD) const {
    if (dump) cout << "calc_f_cp" << endl;
    // Calculate F//
    auto Y0 = xiD * sqrt(m_pars.Kprf * m_pars.Knrf) *
            (m_pars.Crf * m_cache.at("sin(2beta)") -
             m_pars.Srf * m_cache.at("cos(2beta)"));
    if (m_rb < epsilon) return Y0;
    auto Xplus = m_rb * sqrt(m_pars.Kprf * m_pars.Kpwf) *
            (m_pars.Cpp * m_cache.at("sin(2beta + gamma)") -
             m_pars.Spp * m_cache.at("cos(2beta + gamma)"));
    auto Xmins = m_rb * sqrt(m_pars.Knrf * m_pars.Knwf) *
            (-m_pars.Cpn * m_cache.at("sin(2beta + gamma)") -
              m_pars.Spn * m_cache.at("cos(2beta + gamma)"));
    auto Y2 = xiD * m_cache.at("rb2") * sqrt(m_pars.Kpwf * m_pars.Knwf) *
            (m_pars.Cwf * m_cache.at("sin(2beta + 2gamma)") -
             m_pars.Swf * m_cache.at("cos(2beta + 2gamma)"));
    return Y0 + Xplus - Xmins + Y2;
}

DDBPars::ddpair DDBPars::calc_ud_dh() const {
    if (dump) cout << "calc_ud_dh" << endl;
    // Calculate U and D //
    auto u = 0.5 * (1. + m_cache.at("rb2")) * (m_pars.Knrf + m_pars.Kprf);
    auto d = 0.5 * (1. - m_cache.at("rb2")) * (m_pars.Knrf - m_pars.Kprf);
    if (m_rb >= epsilon) {
        const auto sqrtkk = sqrt(m_pars.Kprf * m_pars.Knrf);
        u += 2. * m_rb * sqrtkk * m_cache.at("cos(delb)") *
                (m_pars.C * m_cache.at("cos(gamma)") -
                 m_pars.S * m_cache.at("sin(gamma)"));
        d += 2. * m_rb * sqrtkk * m_cache.at("sin(delb)") *
                (m_pars.S * m_cache.at("cos(gamma)") -
                 m_pars.C * m_cache.at("sin(gamma)"));
    }
    return make_pair(u, d);
}

double DDBPars::calc_f_dh() const {
    if (dump) cout << "calc_f_dh" << endl;
    // Calculate F//
    auto f0 = sqrt(m_pars.Kprf * m_pars.Knrf) *
            (m_pars.C * m_cache.at("sin(2beta)") +
             m_pars.S * m_cache.at("cos(2beta)"));
    if (m_rb < epsilon) return f0;
    auto prod1 = m_cache.at("sin(2beta + gamma)") * m_cache.at("cos(delb)");
    auto prod2 = m_cache.at("cos(2beta + gamma)") * m_cache.at("sin(delb)");
    auto f11 = m_pars.Knrf * (prod1 + prod2);
    auto f12 = m_pars.Kprf * (prod1 - prod2);
    auto f2 = -m_cache.at("rb2") *
            (m_pars.S * m_cache.at("cos(2beta + 2gamma)") -
             m_pars.C * m_cache.at("sin(2beta + 2gamma)"));
    return f0 + m_rb * (f11 + f12) + f2;
}

DDBPars::ddpair DDBPars::calc_ud_dhcp(int xiD) const {
    if (dump) cout << "calc_ud_dhcp" << endl;
    if (m_rb < epsilon) return make_pair(1., 0.);
    // Calculate U and D//
    auto u = 1. + m_cache.at("rb2") + 2. * xiD * m_rb *
             m_cache.at("cos(delb)") * m_cache.at("cos(gamma)");
    auto d = -2. * xiD * m_rb *
             m_cache.at("sin(delb)") * m_cache.at("sin(gamma)");
    return make_pair(u, d);
}

double DDBPars::calc_f_dhcp(int xiD) const {
    if (dump) cout << "calc_f_dhcp" << endl;
    // Calculate F//
    auto f0 = xiD * m_cache.at("sin(2beta)");
    if (m_rb < epsilon) return f0;
    auto f1 = -2. * m_rb * m_cache.at("cos(delb)") * m_cache.at("sin(2beta + gamma)");
    auto f2 = m_cache.at("rb2") * xiD * m_cache.at("sin(2beta + 2gamma)");
    return f0 + f1 + f2;
}

DDBPars::ddpair DDBPars::coefs_dd() const {
    auto ud = calc_ud();
    auto f = calc_f();
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDBPars::ddpair DDBPars::coefs_cp(int xiD) const {
    auto ud = calc_ud_cp(xiD);
    auto f = calc_f_cp(xiD);
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDBPars::ddpair DDBPars::coefs_dh() const {
    auto ud = calc_ud_dh();
    auto f = calc_f_dh();
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDBPars::ddpair DDBPars::coefs_dhcp(int xiD) const {
    auto ud = calc_ud_dhcp(xiD);
    auto f = calc_f_dhcp(xiD);
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDBPars::ddpair DDBPars::coefs(int16_t bbin, int16_t dbin, dtypes dt) const {
    set_pars(bbin, dbin);
    switch (dt) {
    case dtypes::KsPIPI: return coefs_dd();
    case dtypes::CPn:    return coefs_cp(-1);
    case dtypes::CPp:    return coefs_cp(1);
    case dtypes::Dh:     return coefs_dh();
    case dtypes::DhCPn:  return coefs_dhcp(-1);
    case dtypes::DhCPp:  return coefs_dhcp(1);
    default: return make_pair(1., 0.);
    }
}

DDBPars::ddpair DDBPars::rawud(int16_t bbin, int16_t dbin, dtypes dt) const {
    if (dump) cout << "rawud" << endl;
    set_pars(bbin, dbin);
    switch (dt) {
    case dtypes::KsPIPI: return calc_ud();
    case dtypes::CPn:    return calc_ud_cp(-1);
    case dtypes::CPp:    return calc_ud_cp(1);
    case dtypes::Dh:     return calc_ud_dh();
    case dtypes::DhCPn:  return calc_ud_dhcp(-1);
    case dtypes::DhCPp:  return calc_ud_dhcp(1);
    default: return make_pair(1., 0.);
    }
}

DDBPars::Params::Params() :
    Kp(1.), Kn(1.), C(0.), S(0.),
    Kprf(1.), Knrf(1.), Kpwf(1.), Knwf(1.),
    Crf(0.), Srf(0.), Cwf(0.), Swf(0.),
    Ctp(0.), Stp(0.), Ctn(0.), Stn(0.),
    Cpp(0.), Spp(0.), Cpn(0.), Spn(0.) {}
