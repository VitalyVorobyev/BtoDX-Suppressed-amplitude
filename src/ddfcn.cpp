#include "ddfcn.h"

#include <iostream>
#include <map>
#include <cmath>

#include "cfg.h"

using std::cout;
using std::endl;
using std::log;

using std::vector;
using vectd = vector<double>;
using Events = vector<DBEvt>;
using DDMap = std::unordered_map<int16_t,
              std::unordered_map<int16_t, std::pair<double, double>>>;

DDFcn::DDFcn(fitmode fmode, const Events& evtv, uint16_t sb, AbsDDPars& pars) :
   m_mode(fmode), m_evtv(evtv), m_single_bin(sb), m_pars(pars) {}

double DDFcn::operator()(const vectd& par) const {
    m_pars.set_beta(par[0]);

    auto pdf = Cfg::pdf();
    DDMap csmap;  // cached coefficients
    for (auto bbin = -8; bbin <= 8; bbin++) if (bbin)
        for (auto dbin = -8; dbin <= 8; dbin++) if (dbin)
            csmap[bbin][dbin] = m_pars.coefs(bbin, dbin);

    uint64_t cnt = 0;
    uint64_t badcnt = 0;
    double loglh = 0.;
    for (const auto& evt : m_evtv) {
        if (!m_single_bin || (m_single_bin == abs(evt.Bbin()))) {
            auto& cs = csmap[evt.Bbin()][evt.Dbin()];
            auto pdfv = pdf(evt.t(), evt.Tag(), cs.first, cs.second);
            if (pdfv > 0) {
                loglh += log(pdfv);
                cnt++;
            } else {
                loglh -= 100 * pdfv;
                badcnt++;
            }
        }
    }
    loglh *= -2.;
    cout << "llh: " << loglh
         << ", beta: " << par[0]
         << ", cnt: " << cnt
         << ", badcnt: " << badcnt
         << endl;
    return loglh;
}
