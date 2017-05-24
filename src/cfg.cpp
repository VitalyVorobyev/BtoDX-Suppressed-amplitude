#include "cfg.h"

#include <iostream>
#include <cmath>

using std::endl;
using std::cout;

double Cfg::mean = 0.;
double Cfg::sigma = 0.;
double Cfg::fbkg = 0.;
double Cfg::wtag = 0.;

int Cfg::ncp = 4 * std::pow(10, 4);
int Cfg::nfl = 5 * std::pow(10, 5);

std::string Cfg::data_path("/home/vitaly/B0toD0pipi/btoucbard/data/");
std::string Cfg::data_file("data.txt");
std::string Cfg::posi_cp_file("data_posi_cp.txt");
std::string Cfg::nega_cp_file("data_nega_cp.txt");
std::string Cfg::kpi_file("data_kpi.txt");
std::string Cfg::pik_file("data_pik.txt");

// CPV and b -> c ubar d parameters
double Cfg::rd_cp = 1.;
double Cfg::rd_kpi = 0.063;
double Cfg::deld_nega_cp = 180.;  // deg
double Cfg::deld_posi_cp = 0.;  // deg
double Cfg::deld_kpi = 10.;  // deg
double Cfg::rb = 0.05;
double Cfg::delb = 0.;  // deg
double Cfg::beta = 23.;  // deg
double Cfg::ckmgamma = 71.;  // deg

double Cfg::rd_pik() {
    return rd_kpi > 0 ? 1. / rd_kpi : 10000;
}

void Cfg::print_config() {
    cout << "rd " << rd_kpi << endl
         << "deld " << deld_kpi << endl
         << "rb " << rb << endl
         << "delb " << delb << endl
         << "beta " << beta << endl
         << "gamma " << ckmgamma << endl
         << "mean " << mean << endl
         << "sigma " << sigma << endl
         << "fbkg " << fbkg << endl
         << "wtag " << wtag << endl;
}

libTatami::ToyPdf Cfg::pdf() {
    return libTatami::ToyPdf(mean, sigma, fbkg, wtag);
}

int Cfg::n_posi_cp() {return ncp / 2;}
int Cfg::n_nega_cp() {return ncp / 2;}
int Cfg::n_flv_rs() {return nfl / 2;}
int Cfg::n_flv_ws() {return nfl * rd_kpi / 2;}

Lambda_b2cud Cfg::lamf(dtypes type) {
    switch (type) {
    case dtypes::CPp:
        return Lambda_b2cud(rd_cp, deld_posi_cp, rb, delb, beta, ckmgamma);
    case dtypes::CPn:
        return Lambda_b2cud(rd_cp, deld_nega_cp, rb, delb, beta, ckmgamma);
    case dtypes::KPI:
        return Lambda_b2cud(rd_kpi, deld_kpi, rb, delb, beta, ckmgamma);
    case dtypes::PIK:
        return Lambda_b2cud(rd_pik(), deld_kpi, rb, delb, beta, ckmgamma);
    default:
        cout << "Cfg::lamf: Unknown eventtype " << type << endl;
        return Lambda_b2cud(1, 0, 0, 0, beta, ckmgamma);
    }
}

Lambda_b2cud Cfg::lamf0(dtypes type) {
    switch (type) {
    case dtypes::CPp:
        return Lambda_b2cud(rd_cp, deld_posi_cp, 0, delb, beta, ckmgamma);
    case dtypes::CPn:
        return Lambda_b2cud(rd_cp, deld_nega_cp, 0, delb, beta, ckmgamma);
    case dtypes::KPI:
        return Lambda_b2cud(rd_kpi, deld_kpi, 0, delb, beta, ckmgamma);
    case dtypes::PIK:
        return Lambda_b2cud(rd_pik(), deld_kpi, 0, delb, beta, ckmgamma);
    default:
        cout << "Cfg::lamf0: Unknown eventtype " << type << endl;
        return Lambda_b2cud(1, 0, 0, 0, beta, ckmgamma);
    }
}

std::pair<int, int> Cfg::nevts(dtypes type) {
    switch (type) {
    case dtypes::CPp:
        return std::make_pair(n_posi_cp() / 2, n_posi_cp() / 2);
    case dtypes::CPn:
        return std::make_pair(n_nega_cp() / 2, n_nega_cp() / 2);
    case dtypes::KPI:
        return std::make_pair(n_flv_rs(), n_flv_ws());
    case dtypes::PIK:
        return std::make_pair(n_flv_ws(), n_flv_rs());
    default:
        cout << "Cfg::nevts: Unknown eventtype " << type << endl;
        return std::make_pair(0, 0);
    }
}

std::string Cfg::dfile(dtypes type) {
    switch (type) {
    case dtypes::CPp:
        return data_path + posi_cp_file;
    case dtypes::CPn:
        return data_path + nega_cp_file;
    case dtypes::KPI:
        return data_path + kpi_file;
    case dtypes::PIK:
        return data_path + pik_file;
    default:
        cout << "Cfg::dfile: Unknown eventtype " << type << endl;
        return data_path + data_file;
    }
}
