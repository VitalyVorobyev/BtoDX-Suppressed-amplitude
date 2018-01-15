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
                 const str& dcfg, const str& bcfg) : AbsDDPars(dcfg, bcfg, beta) {
    setNewParam("rb", rb);
    setNewParam("gamma", gamma);
    setNewParam("delb", delb);

    cout << "New DDBPars instance: rb = " << rb
         << ", beta = " << beta
         << ", gamma = " << gamma
         << ", delb = " << delb
         << endl;

    setBeta(beta);
    setGamma(gamma);
    setRb(rb);
    setDeltaB(delb);
    updateCache();

    print();
}

void DDBPars::setParam(const std::string& name, double val) {
    AbsDDPars::setParam(name, val);
    if (name == "beta") setBeta(val);
    else if (name == "gamma") setGamma(val);
    else if (name == "rb") setRb(val);
    else if (name == "delb") setDeltaB(val);
}

void DDBPars::setRb(double x) {
    m_cache["rb2"] = x * x;
}

void DDBPars::setDeltaB(double x) {
    m_cache["sin(delb)"] = sin(x);
    m_cache["cos(delb)"] = cos(x);
}

void DDBPars::setBeta(double x) {
    m_cache["cos(2beta)"] = cos(2. * x);
    m_cache["sin(2beta)"] = sin(2. * x);
    updateCache();
}

void DDBPars::setGamma(double x) {
    m_cache["cos(gamma)"] = cos(x);
    m_cache["sin(gamma)"] = sin(x);
    updateCache();
}

void DDBPars::updateCache() {
    const auto& beta = getParam("beta");
    const auto& gamma = getParam("gamma");
    m_cache["cos(2beta + gamma)"] = cos(2. * beta + gamma);
    m_cache["sin(2beta + gamma)"] = sin(2. * beta + gamma);
    m_cache["cos(2beta + 2gamma)"] = cos(2. * (beta + gamma));
    m_cache["sin(2beta + 2gamma)"] = sin(2. * (beta + gamma));
}

void DDBPars::set_pars(int16_t bbin, int16_t dbin) const {
    if (!bbin && !dbin) return;
    if (dump) cout << "Get pars" << endl;
    // D meson binned Dalitz params
    if (dbin > 0) {
        dbin--;  // bin to index
        m_pars.Kp = getInt("K+", dbin);
        m_pars.Kn = getInt("K-", dbin);
        m_pars.C = getInt("C", dbin);
        m_pars.S = getInt("S", dbin);
    } else if (dbin < 0) {
        dbin = -dbin - 1;  // bin to index
        m_pars.Kp = getInt("K-", dbin);
        m_pars.Kn = getInt("K+", dbin);
        m_pars.C = getInt("C", dbin);
        m_pars.S = -getInt("S", dbin);
    }
    if (bbin > 0) {
        bbin--;
        m_pars.Kprf = getInt("K+rf", bbin);
        m_pars.Knrf = getInt("K-rf", bbin);
        m_pars.Crf = getInt("Crf", bbin);
        m_pars.Srf = getInt("Srf", bbin);
    } else if (bbin < 0) {
        bbin = -bbin - 1;
        m_pars.Kprf = getInt("K-rf", bbin);
        m_pars.Knrf = getInt("K+rf", bbin);
        m_pars.Crf = getInt("Crf", bbin);
        m_pars.Srf = -getInt("Srf", bbin);
    }
    const auto rb = getParam("rb");
    if (rb < epsilon) return;
    // B meson binned Dalitz params
    if (bbin > 0) {
        m_pars.Kpwf = getInt("K+wf", bbin);
        m_pars.Knwf = getInt("K-wf", bbin);
        m_pars.Cwf = getInt("Cwf", bbin);
        m_pars.Swf = getInt("Swf", bbin);

        m_pars.Ctp = getInt("Ctp", bbin);
        m_pars.Stp = getInt("Stp", bbin);
        m_pars.Ctn = getInt("Ctn", bbin);
        m_pars.Stn = getInt("Stn", bbin);

        m_pars.Cpp = getInt("Cpp", bbin);
        m_pars.Spp = getInt("Spp", bbin);
        m_pars.Cpn = getInt("Cpn", bbin);
        m_pars.Spn = getInt("Spn", bbin);
    } else if (bbin < 0) {
        m_pars.Kpwf = getInt("K-wf", bbin);
        m_pars.Knwf = getInt("K+wf", bbin);
        m_pars.Cwf = getInt("Cwf", bbin);
        m_pars.Swf = -getInt("Swf", bbin);

        m_pars.Ctp = getInt("Ctn", bbin);
        m_pars.Stp = getInt("Stn", bbin);
        m_pars.Ctn = getInt("Ctp", bbin);
        m_pars.Stn = getInt("Stp", bbin);

        m_pars.Cpp = getInt("Cpn", bbin);
        m_pars.Spp = getInt("Spn", bbin);
        m_pars.Cpn = getInt("Cpp", bbin);
        m_pars.Spn = getInt("Spp", bbin);
    }
}

DDBPars::ddpair DDBPars::calc_ud() const {
    if (dump) cout << "calc_ud" << endl;
    // Calculate U //
    auto KpKn = m_pars.Kp * m_pars.Kn;
    auto p11 = m_pars.Kn * m_pars.Kprf;
    auto p12 = m_pars.Kp * m_pars.Knrf;
    auto p1 = 0.5 * (p11 + p12);

    const auto rb = getParam("rb");
    if (rb < epsilon)
        return make_pair(p1, 0.5 * (p11 - p12));

    auto p20 = rb * sqrt(KpKn * m_pars.Kprf * m_pars.Kpwf);
    auto p21 = (m_pars.Ctp * m_pars.C +
                m_pars.Stp * m_pars.S) * cache("cos(gamma)");
    auto p22 = (m_pars.Ctp * m_pars.S +
                m_pars.Stp * m_pars.C) * cache("sin(gamma)");
    auto p2 = p20 * (p21 - p22);

    auto p30 = rb * sqrt(KpKn * m_pars.Knrf * m_pars.Knwf);
    auto p31 = (m_pars.Ctn * m_pars.C -
                m_pars.Stn * m_pars.S) * cache("cos(gamma)");
    auto p32 = (m_pars.Ctn * m_pars.S +
                m_pars.Stn * m_pars.C) * cache("sin(gamma)");
    auto p3 = p30 * (p31 - p32);

    auto p41 = m_pars.Kp * m_pars.Kpwf;
    auto p42 = m_pars.Kn * m_pars.Knwf;
    auto p4 = 0.5 * cache("rb2") * (p41 + p42);
    auto u = p1 + p2 + p3 + p4;

    // Calcuate D //
    p1 = 0.5 * (p11 - p12);
    p4 = 0.5 * cache("rb2") * (p41 - p42);
    auto d = p1 + p2 - p3 + p4;
    return make_pair(u, d);
}

double DDBPars::calc_f() const {
    if (dump) cout << "calc_f" << endl;
    // Calculate F //
    auto KpKn = m_pars.Kp * m_pars.Kn;
    auto p10 = sqrt(KpKn * m_pars.Kprf * m_pars.Knrf);
    auto p11 = (m_pars.Srf * m_pars.C -
                m_pars.Crf * m_pars.S) * cache("cos(2beta)");
    auto p12 = (m_pars.Crf * m_pars.C +
                m_pars.Srf * m_pars.S) * cache("sin(2beta)");
    auto p1 = -p10 * (p11 - p12);

    const auto rb = getParam("rb");
    if (rb < epsilon) return p1;

    auto p20 = rb * m_pars.Kn * sqrt(m_pars.Kprf * m_pars.Knwf);
    auto p21 = m_pars.Cpp * cache("sin(2beta + gamma)");
    auto p22 = m_pars.Spp * cache("cos(2beta + gamma)");
    auto p2 = p20 * (p21 - p22);

    auto p30 = rb * m_pars.Kp * sqrt(m_pars.Knrf * m_pars.Kpwf);
    auto p31 = m_pars.Cpn * cache("sin(2beta + gamma)");
    auto p32 = m_pars.Spn * cache("cos(2beta + gamma)");
    auto p3 = p30 * (p31 + p32);

    auto p40 = 0.5 * m_cache.at("rb2") *
            sqrt(KpKn * m_pars.Kpwf * m_pars.Knwf);
    auto p41 = (m_pars.Swf * m_pars.C +
                m_pars.Cwf * m_pars.S) * cache("cos(2beta + 2gamma)");
    auto p42 = (m_pars.Cwf * m_pars.C -
                m_pars.Swf * m_pars.S) * cache("sin(2beta + 2gamma)");
    auto p4 = -p40 * (p41 - p42);
    return p1 + p2 + p3 + p4;
}

DDBPars::ddpair DDBPars::calc_ud_cp(int xiD) const {
    if (dump) cout << "calc_ud_cp" << endl;
    // Calculate U and D //
    auto Xplus = 0.5 * (m_pars.Kprf + cache("rb2") * m_pars.Kpwf);
    auto Xmins = 0.5 * (m_pars.Knrf + cache("rb2") * m_pars.Knwf);

    const auto rb = getParam("rb");
    if (rb >= epsilon) {
        Xplus += xiD * rb * sqrt(m_pars.Kprf * m_pars.Kpwf) *
                (m_pars.Ctp * cache("cos(gamma)") -
                 m_pars.Stp * cache("sin(gamma)"));
        Xmins += xiD * rb * sqrt(m_pars.Knrf * m_pars.Knwf) *
                (m_pars.Ctn * cache("cos(gamma)") +
                 m_pars.Stn * cache("sin(gamma)"));
    }
    return make_pair(Xplus + Xmins, Xplus - Xmins);
}

double DDBPars::calc_f_cp(int xiD) const {
    if (dump) cout << "calc_f_cp" << endl;
    // Calculate F//
    auto Y0 = xiD * sqrt(m_pars.Kprf * m_pars.Knrf) *
            (m_pars.Crf * cache("sin(2beta)") -
             m_pars.Srf * cache("cos(2beta)"));

    const auto rb = getParam("rb");
    if (rb < epsilon) return Y0;
    auto Xplus = rb * sqrt(m_pars.Kprf * m_pars.Kpwf) *
            (m_pars.Cpp * cache("sin(2beta + gamma)") -
             m_pars.Spp * cache("cos(2beta + gamma)"));
    auto Xmins = rb * sqrt(m_pars.Knrf * m_pars.Knwf) *
            (-m_pars.Cpn * cache("sin(2beta + gamma)") -
              m_pars.Spn * cache("cos(2beta + gamma)"));
    auto Y2 = xiD * cache("rb2") * sqrt(m_pars.Kpwf * m_pars.Knwf) *
            (m_pars.Cwf * cache("sin(2beta + 2gamma)") -
             m_pars.Swf * cache("cos(2beta + 2gamma)"));
    return Y0 + Xplus - Xmins + Y2;
}

DDBPars::ddpair DDBPars::calc_ud_dh() const {
    if (dump) cout << "calc_ud_dh" << endl;
    // Calculate U and D //
    auto u = 0.5 * (1. + cache("rb2")) * (m_pars.Kn + m_pars.Kp);
    auto d = 0.5 * (1. - cache("rb2")) * (m_pars.Kn - m_pars.Kp);

    if (dump) cout << "Zero rb: " << endl
                   << "  u = " << u << endl
                   << "  d = " << d << endl;

    const auto rb = getParam("rb");
    if (rb >= epsilon) {
        const auto sqrtkk = sqrt(m_pars.Kp * m_pars.Kn);
        u += 2. * rb * sqrtkk * cache("cos(delb)") *
                (m_pars.C * cache("cos(gamma)") +
                 m_pars.S * cache("sin(gamma)"));
        d += 2. * rb * sqrtkk * m_cache.at("sin(delb)") *
                (m_pars.S * cache("cos(gamma)") -
                 m_pars.C * cache("sin(gamma)"));
    }

    if (dump) cout << "Non-zero rb: " << endl
                   << "  u = " << u << endl
                   << "  d = " << d << endl;
    return make_pair(u, d);
}

double DDBPars::calc_f_dh() const {
    if (dump) cout << "calc_f_dh" << endl;
    // Calculate F//
    auto f0 = sqrt(m_pars.Kp * m_pars.Kn) *
            (m_pars.C * cache("sin(2beta)") +
             m_pars.S * cache("cos(2beta)"));
    const auto rb = getParam("rb");
    if (rb < epsilon) return f0;
    auto prod1 = cache("sin(2beta + gamma)") * cache("cos(delb)");
    auto prod2 = cache("cos(2beta + gamma)") * cache("sin(delb)");
    auto f11 = m_pars.Kp * (prod1 + prod2);
    auto f12 = m_pars.Kn * (prod1 - prod2);
    auto f2 = cache("rb2") *
            (m_pars.C * cache("sin(2beta + 2gamma)") -
             m_pars.S * cache("cos(2beta + 2gamma)"));
    if (dump) cout << " f0:  " << f0 << endl
                   << " f11: " << f11 << endl
                   << " f12: " << f12 << endl
                   << " f1:  " << rb * (f11 + f12) << endl
                   << " f2:  " << f2 << endl
                   << " rb:  " << rb << endl
                   << " rb2: " << cache("rb2") << endl
                   << " f : " << f0 + rb * (f11 + f12) + f2 << endl;
    if (dump) cout << "f2 = " << cache("rb2") << " * (" << m_pars.C << " * "
                   << cache("sin(2beta + 2gamma)") << " - " << m_pars.S << " * "
                   << cache("cos(2beta + 2gamma)") << ")" << endl;
    return f0 + rb * (f11 + f12) + f2;
}

DDBPars::ddpair DDBPars::calc_ud_dhcp(int xiD) const {
    if (dump) cout << "calc_ud_dhcp" << endl;
    const auto rb = getParam("rb");
    if (rb < epsilon) return make_pair(1., 0.);
    // Calculate U and D//
    auto u = 1. + cache("rb2") + 2. * xiD * rb *
             cache("cos(delb)") * cache("cos(gamma)");
    auto d = -2. * xiD * rb *
             cache("sin(delb)") * cache("sin(gamma)");
    return make_pair(u, d);
}

double DDBPars::calc_f_dhcp(int xiD) const {
    if (dump) cout << "calc_f_dhcp" << endl;
    // Calculate F//
    auto f0 = xiD * cache("sin(2beta)");

    const auto rb = getParam("rb");
    if (rb < epsilon) return f0;
    auto f1 = 2. * rb * cache("cos(delb)") * cache("sin(2beta + gamma)");
    auto f2 = cache("rb2") * xiD * cache("sin(2beta + 2gamma)");
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
    auto c = ud.second / ud.first;
    auto s = f / ud.first;
    if (std::fabs(s) > 1.) {
        static int flag = 0;
        if (flag++ < 16)
            cerr << "WARNING: |s| > 1" << endl;
        s = s > 0 ? 1. : -1.;
    }
    return make_pair(c, s);
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
    case dtypes::DKs:    return coefs_dh();
    case dtypes::DCPnKs: return coefs_dhcp(-1);
    case dtypes::DCPpKs: return coefs_dhcp(1);
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
    case dtypes::DKs:    return calc_ud_dh();
    case dtypes::DCPnKs: return calc_ud_dhcp(-1);
    case dtypes::DCPpKs: return calc_ud_dhcp(1);
    default: return make_pair(1., 0.);
    }
}

DDBPars::Params::Params() :
    Kp(1.), Kn(1.), C(0.), S(0.),
    Kprf(1.), Knrf(1.), Kpwf(1.), Knwf(1.),
    Crf(0.), Srf(0.), Cwf(0.), Swf(0.),
    Ctp(0.), Stp(0.), Ctn(0.), Stn(0.),
    Cpp(0.), Spp(0.), Cpn(0.), Spn(0.) {}
