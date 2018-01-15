#pragma once

#include <memory>
#include <cstdint>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <set>

#include "fitmodes.h"
#include "cfg.h"

// Forward declarations
class AbsDDPars;
class BEvt;
class DBEvt;
namespace libTatami {
    class ToyPdf;
}

namespace ROOT {
    namespace Math {
        class Minimizer;
    }
}

class CPVMinimizer {
    using MapCS = std::map<int16_t,
                  std::map<int16_t,
                  std::pair<double, double>>>;
    /** Binnned Dalitz C and S keeper */
    using SMapCS = std::map<int16_t,
                   std::map<dtypes,
                   std::pair<double, double>>>;
    using ParsLookUp = std::map<std::string, uint16_t>;

    // Static consts
    static constexpr double rad_to_deg = 180. / M_PI;
    static const std::map<fitmode, bool> fitphase;
    static const uint16_t m_nbins = 8;
    static ParsLookUp m_parsLookUp;
    static std::set<std::string> m_fixed_vars;

    static std::string m_min_name;
    static std::string m_min_alg;
    static uint16_t m_single_bin;
    static bool m_rndm_init;
    static fitmode m_mode;
    static std::unique_ptr<libTatami::ToyPdf> m_pdf;
    static bool m_dh_flag;

    static std::unique_ptr<AbsDDPars> m_pars;
    /** B0 -> Dcp pi+ pi- events */
    static std::vector<BEvt> sevtv;
    /** B0 -> (Ks0 pi+ pi-)_D0 pi+ pi- events */
    static std::vector<DBEvt> sdevtv;

    static double lhfcn(const double* x);
    static double DDLH();
    static double CPLH();
    static double DhLH();

    static uint16_t addPar(ROOT::Math::Minimizer& m,
                           const std::string& name, double init, double var);
    static uint16_t m_parIdx;
    static uint16_t fixPars(ROOT::Math::Minimizer& m);

    static Cfg::ExpSetup m_expCfg;

    // private constructor
    CPVMinimizer() {}

 public:
    static void MakeFit(std::unique_ptr<AbsDDPars>& pars,
                        fitmode mode, uint16_t sb=0);

    static void fixPar(const std::string& name);
    // Data management
    static void FlushData();
    static void ReadData(const std::string &input_file_name);
    static void ReadBData(const std::string &input_file_name);
    static void setSetup(Cfg::ExpSetup setup) {
        m_expCfg = setup;
    }
};
