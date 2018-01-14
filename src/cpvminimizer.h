#pragma once

#include <memory>
#include <cstdint>
#include <string>
#include <map>
#include <cmath>
#include <vector>

#include "fitmodes.h"

// Forward declarations
class AbsDDPars;
class BEvt;
class DBEvt;
namespace libTatami {
    class ToyPdf;
}

class CPVMinimizer {
    using MapCS = std::map<int16_t,
                  std::map<int16_t,
                  std::pair<double, double>>>;
    /** Binnned Dalitz C and S keeper */
    using SMapCS = std::map<int16_t,
                   std::map<dtypes,
                   std::pair<double, double>>>;
    // Static consts
    static constexpr double rad_to_deg = 180. / M_PI;
    static const std::map<fitmode, bool> fitphase;
    static const uint16_t m_nbins = 8;

    static std::string m_min_name;
    static std::string m_min_alg;
    static uint16_t m_single_bin;
    static bool m_rndm_init;
    static fitmode m_mode;
    static libTatami::ToyPdf m_pdf;
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

    // private constructor
    CPVMinimizer() {}

 public:
    static void MakeFit(std::unique_ptr<AbsDDPars>& pars,
                        fitmode mode, uint16_t sb=0);
    // Data management
    static void FlushData();
    static void ReadData(const std::string &input_file_name);
    static void ReadBData(const std::string &input_file_name);
};
