#pragma once

#include <string>
#include <utility>  // pair
#include <memory>
#include <map>

#include "fitmodes.h"

// Forward declarations
class DDBPars;
class DDMPars;
namespace libTatami {
    class ToyPdf;
}

/** @brief The Cfg class keeps configuration */
class Cfg {
 public:
   enum class ExpSetup : int {Belle = 0, LHCb = 1, Perfect = 2};

 private:
    Cfg();

    class ExpCfg {
        double m_mean;
        double m_sigma;
        double m_wtag;
        double m_fbkg;

     public:
        ExpCfg(double m, double s, double w, double f) :
            m_mean(m), m_sigma(s), m_wtag(w), m_fbkg(f) {}

        double mean() const {return m_mean;}
        double sigma() const {return m_sigma;}
        double wtag() const {return m_wtag;}
        double fbkg() const {return m_fbkg;}
    };

    const static ExpCfg belleCfg;
    const static ExpCfg lhcbCfg;
    const static ExpCfg perfCfg;

    const static std::map<ExpSetup, ExpCfg> expCfgMap;

    // Experimental conditions
    static int ncp;
    static int nfl;
    static inline int n_posi_cp();
    static inline int n_nega_cp();
    static inline int n_flv_rs();
    static inline int n_flv_ws();

    // CPV and b -> c ubar d parameters
    static double rd_cp;
    static double rd_kpi;
    static inline double rd_pik();
    static double deld_nega_cp;  // deg
    static double deld_posi_cp;  // deg
    static double deld_kpi;  // deg
    static double rb;
    static double delb;  // deg
    static double beta;  // deg
    static double ckmgamma;  // deg
    static double charm_x;
    static double charm_y;
    static double dtlim;

 public:
    static void print_config();

    static std::pair<int, int> nevts(dtypes type);
    static std::string dfile(dtypes type, uint32_t nevt=0);
    static std::unique_ptr<libTatami::ToyPdf> pdf(ExpSetup exp);

    static void set_beta(double x) {beta = x;}
    static void set_dtlim(double x) {dtlim = x;}
    static void set_rb(double x) {rb = x;}
    static void set_delb(double x) {delb = x;}
    static void set_ncp(double x) {ncp = x;}
    static void set_nfl(double x) {nfl =x;}
    static void set_charm_mix(double x, double y) {
        charm_x = x; charm_y = y;
    }
    static double get_rb() {return rb;}
    static double get_x() {return charm_x;}
    static double get_y() {return charm_y;}

    static std::unique_ptr<DDBPars> wfpars(const std::string& dcfg,
                                           const std::string& bcfg);
    static std::unique_ptr<DDBPars> wfpars(uint16_t seed, uint16_t idx);
    static std::unique_ptr<DDMPars> cmpars(const std::string& dcfg,
                                           const std::string& bcfg);
    static std::unique_ptr<DDMPars> cmpars();
    /** B -> D0 pi+ pi- config file name */
    static std::string get_bcfg(uint16_t seed, uint16_t idx);
    /** D -> Ks0 pi+ pi- config file name */
    static std::string get_dcfg();
    /** Data file name for rB != 0 events */
    static std::string wfdtdist(uint16_t seed, uint16_t idx,
                                double rb, dtypes type);
    /** Data file name for events with charm mixing */
    static std::string cmdtdist(double x, double y, dtypes type);
};
