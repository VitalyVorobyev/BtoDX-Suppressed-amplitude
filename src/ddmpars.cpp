#include "ddmpars.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using std::string;
using std::vector;
using std::sin;
using std::cos;

using std::cout;
using std::cerr;
using std::endl;
using std::ifstream;

using std::pair;
using std::make_pair;

constexpr bool dump = false;
constexpr double epsilon = 0.0001;

DDMPars::DDMPars(double beta, double x, double y,
                 const string& dcfg, const string& bcfg) :
    AbsDDPars(dcfg, bcfg, beta) {
    setNewParam("x", x);
    setNewParam("y", y);
    setBeta(beta);
}

bool DDMPars::zeroMixing() const {
    return (std::fabs(getParam("x")) < epsilon) &&
           (std::fabs(getParam("y")) < epsilon);
}

void DDMPars::setParam(const std::string &name, double val) {
    AbsDDPars::setParam(name, val);
    if (name == "beta") setBeta(val);
}

void DDMPars::setBeta(double x) {
    m_cache["sin(2beta)"] = sin(2. * x);
    m_cache["cos(2beta)"] = cos(2. * x);
}

void DDMPars::set_pars(int16_t bbin, int16_t dbin) const {
    if (dump) cout << "Set pars" << endl;
    if (!bbin && !dbin) return;
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
    } else if (bbin < 0){
        bbin = -bbin - 1;
        m_pars.Kprf = getInt("K-rf", bbin);
        m_pars.Knrf = getInt("K+rf", bbin);
        m_pars.Crf = getInt("Crf", bbin);
        m_pars.Srf = -getInt("Srf", bbin);
    }
}

DDMPars::ddpair DDMPars::calc_ud() const {
    const auto x = getParam("x");
    const auto y = getParam("y");
    auto srhKpKn = 0.5 * sqrt(m_pars.Kp * m_pars.Kn);
    auto p1 = 0.5 * m_pars.Kprf * m_pars.Kn;
    auto p2 = 0.5 * m_pars.Knrf * m_pars.Kp;

    if (zeroMixing()) return make_pair(p1 + p2, p1 - p2);

    auto p3 = (y * m_pars.C + x * m_pars.S) * m_pars.Kprf;
    auto p4 = (y * m_pars.C - x * m_pars.S) * m_pars.Knrf;

    auto u = p1 + p2 + srhKpKn * (p3 + p4);
    auto d = p1 - p2 + srhKpKn * (p3 - p4);
    return make_pair(u, d);
}

DDMPars::ddpair DDMPars::calc_ud_dcp(int16_t xid) const {
    if (dump) cout << "DDMPars::calc_ud_dcp" << endl;
    return make_pair(
                0.5 * (m_pars.Kprf + m_pars.Knrf) * (1. + xid * getParam("y")),
                0.5 * (m_pars.Kprf - m_pars.Knrf) * (1. + xid * getParam("y"))
                );
}

DDMPars::ddpair DDMPars::calc_ud_dh() const {
    if (dump) cout << "DDMPars::calc_ud_dh" << endl;
    const auto x = getParam("x");
    const auto y = getParam("y");
    auto srhKpKn = sqrt(m_pars.Kp * m_pars.Kn);
    return make_pair(
                0.5 * (m_pars.Kp + m_pars.Kn) + y * m_pars.C * srhKpKn,
                0.5 * (m_pars.Kp - m_pars.Kn) + x * m_pars.S * srhKpKn);
}

DDMPars::ddpair DDMPars::calc_ud_dh_dcp(int16_t xid) const {
    if (dump) cout << "DDMPars::calc_ud_dh_dcp" << endl;
    return make_pair(1. + xid * getParam("y"), 0.);
}

double DDMPars::calc_f() const {
    const auto x = getParam("x");
    const auto y = getParam("y");
    auto KKrf = m_pars.Kprf * m_pars.Knrf;
    auto p00 = sqrt(KKrf * m_pars.Kp * m_pars.Kn);
    auto p01 = (m_pars.C * m_pars.Crf +
                m_pars.S * m_pars.Srf) * cache("sin(2beta)");
    auto p02 = (m_pars.S * m_pars.Crf -
                m_pars.C * m_pars.Srf) * cache("cos(2beta)");
    auto p0 = p00 * (p01 + p02);

    if (zeroMixing()) return p0;

    auto p11 = y * (m_pars.Srf * cache("cos(2beta)") -
                    m_pars.Crf * cache("sin(2beta)")) *
            (m_pars.Kp + m_pars.Kn);
    auto p12 = x * (m_pars.Crf * cache("cos(2beta)") +
                    m_pars.Srf * cache("sin(2beta)")) *
            (m_pars.Kp - m_pars.Kn);
    auto p1 = +0.5 * sqrt(KKrf) * (p11 + p12);
    return p0 + p1;
}

double DDMPars::calc_f_dcp(int16_t xid) const {
    return xid * sqrt(m_pars.Kprf * m_pars.Knrf) *
            (m_pars.Crf * cache("sin(2beta)") -
             m_pars.Srf * cache("cos(2beta)")) * (1. - xid * getParam("y"));
}

double DDMPars::calc_f_dh(int16_t xih) const {
    if (dump) cout << "DDMPars::calc_f_dh" << endl;
    const auto x = getParam("x");
    const auto y = getParam("y");
    return xih * sqrt(m_pars.Kp * m_pars.Kn) *
            (m_pars.C * cache("sin(2beta)") - m_pars.S * cache("cos(2beta)"))
            + 0.5 * xih * (x * cache("cos(2beta)") * (m_pars.Kp - m_pars.Kn) -
                           y * cache("sin(2beta)") * (m_pars.Kp + m_pars.Kn));
}

double DDMPars::calc_f_dh_dcp(int16_t xid, int16_t xih) const {
    return xih * xid * cache("sin(2beta)") * (1. - xid * getParam("y"));
}

DDMPars::ddpair DDMPars::coefs_dd() const {
    auto ud = calc_ud();
    auto f = calc_f();
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDMPars::ddpair DDMPars::coefs_cp(int16_t xid) const {
    auto ud = calc_ud_dcp(xid);
    auto f = calc_f_dcp(xid);
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDMPars::ddpair DDMPars::coefs_dh() const {
    auto ud = calc_ud_dh();
    auto f = calc_f_dh(1);
    if (dump) cout << "u " << ud.first << ", d " << ud.second
                   << ", f " << f << endl;
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDMPars::ddpair DDMPars::coefs_dhcp(int16_t xid) const {
    auto ud = calc_ud_dh_dcp(xid);
    auto f = calc_f_dh_dcp(xid, 1);
    if (dump) cout << "u " << ud.first << ", d " << ud.second
                   << ", f " << f << endl;
    return make_pair(ud.second / ud.first, f / ud.first);
}

DDMPars::ddpair DDMPars::coefs(int16_t bbin, int16_t dbin, dtypes dt) const {
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

DDMPars::ddpair DDMPars::rawud(int16_t bbin, int16_t dbin, dtypes dt) const {
    set_pars(bbin, dbin);
    switch (dt) {
    case dtypes::KsPIPI: return calc_ud();
    case dtypes::CPn:    return calc_ud_dcp(-1);
    case dtypes::CPp:    return calc_ud_dcp(1);
    case dtypes::Dh:     return calc_ud_dh();
    case dtypes::DhCPn:  return calc_ud_dh_dcp(1);
    case dtypes::DhCPp:  return calc_ud_dh_dcp(-1);
    case dtypes::DKs:    return calc_ud_dh();
    case dtypes::DCPnKs: return calc_ud_dh_dcp(1);
    case dtypes::DCPpKs: return calc_ud_dh_dcp(-1);
    default: return make_pair(1., 0.);
    }
}
