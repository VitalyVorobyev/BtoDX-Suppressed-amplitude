#pragma once

#include <utility>  // pair
#include <vector>
#include <string>
#include <memory>
#include <cstdint>

#include "fitmodes.h"
#include "cfg.h"

// Forward declarations
class DBEvt;
class BEvt;
class AbsDDPars;

/**
 * @brief The DDGev class generates all types of B0 -> D0 {pi+ pi-, h0} events
 */
class DDGev {
    using pDBEvtVec = std::vector<std::unique_ptr<DBEvt>>;
    using pBEvtVec = std::vector<std::unique_ptr<BEvt>>;

    static constexpr uint16_t m_nbins = 8;
    static std::unique_ptr<pBEvtVec>
    CPDalitz(uint64_t nevt, const AbsDDPars& pars, dtypes type);
    static std::unique_ptr<pBEvtVec>
    DhCP(uint64_t nevt, const AbsDDPars& pars, dtypes type);
    static Cfg::ExpSetup m_expCfg;
    /** Forbid to instantiate this class */
    DDGev() {}

 public:
    /** @brief Generate B0 -> Dcp pi+ pi- events */
    static void CPDalitz(uint64_t ncpp, uint64_t ncpn, const AbsDDPars& pars,
                         const std::string& fname);
    /** @brief Generate B0 -> [D0 -> Ks0 pi+ pi-] pi+ pi- events */
    static void DoubleDalitz(uint64_t nevt, const AbsDDPars& pars,
                             const std::string &fname);
    /** @brief Generate B0 -> Dcp {h0, Ks0} events */
    static void DhCP(uint64_t ncpp, uint64_t ncpn, const AbsDDPars& pars,
                     const std::string& fname);
    /** @brief Generate B0 -> [D0 -> Ks0 pi+ pi-] {h0, Ks0} events */
    static void DhDalitz(uint64_t nevt, const AbsDDPars& pars,
                         const std::string& fname);
    static void setSetup(Cfg::ExpSetup setup) {
        m_expCfg = setup;
    }
};
