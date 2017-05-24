#ifndef CFG_H
#define CFG_H

#include <string>
#include <utility>  // pair

#include "mylibs/libTatami/toypdf.h"

#include "lambda.h"
#include "evt.h"

class Cfg {
 public:
    static void print_config();

    static Lambda_b2cud lamf(dtypes type);
    static Lambda_b2cud lamf0(dtypes type);
    static std::pair<int, int> nevts(dtypes type);
    static std::string dfile(dtypes type);
    static libTatami::ToyPdf pdf();

    static void set_delb(double x) {delb = x;}
    static void set_ncp(double x) {ncp = x;}
    static void set_nfl(double x) {nfl =x;}

 private:
    Cfg();

    // Experimental conditions
    static double mean, sigma, fbkg, wtag;

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

    // Data files
    static std::string data_path;
    static std::string data_file;
    static std::string posi_cp_file;
    static std::string nega_cp_file;
    static std::string kpi_file;
    static std::string pik_file;
};

#endif  // CFG_H
