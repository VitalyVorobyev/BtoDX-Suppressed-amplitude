#include "cfg.h"

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "mylibs/libTatami/toypdf.h"

#include "evt.h"
#include "ddbpars.h"
#include "ddmpars.h"

const auto dtr = M_PI / 180.;
using std::make_unique;

using std::ifstream;
using std::string;

using std::vector;

using std::endl;
using std::cout;

using std::unique_ptr;
using std::to_string;

const Cfg::ExpCfg Cfg::belleCfg(0., 1.25, 0.3, 0.3);
const Cfg::ExpCfg Cfg::lhcbCfg(0., 0.06, 0.08, 0.05);
const Cfg::ExpCfg Cfg::perfCfg(0., 0., 0., 0.);

const std::map<Cfg::ExpSetup, Cfg::ExpCfg> Cfg::expCfgMap {
    {Cfg::ExpSetup::Belle, Cfg::belleCfg},
    {Cfg::ExpSetup::LHCb, Cfg::lhcbCfg},
    {Cfg::ExpSetup::Perfect, Cfg::perfCfg},
};

int Cfg::ncp = 4 * std::pow(10, 4);
int Cfg::nfl = 5 * std::pow(10, 5);

const string pars_path("/home/vitaly/B0toD0pipi/B0toD0pipiFeas/params/");
const string dpars("kspipi_meas_params.txt");

const string data_path("/home/vitaly/B0toD0pipi/btoucbard/data/");
const string data_file("data");
const string posi_cp_file("data_posi_cp");
const string nega_cp_file("data_nega_cp");
const string kpi_file("data_kpi");
const string pik_file("data_pik");
const string kspp_file("data_kspp");

const string kspp_csk_file("csk_tblSymABAC_kspp_1M.txt");
const string kspp_adds_file("adds_tblSymABAC_kspp_1M.txt");
const string kspp_approx_adds_file("adds_approx_tblSymABAC_kspp_1M.txt");

// String labels for data types
const std::map<dtypes, std::string> dtmap = {
    {dtypes::CPp, "cp"},
    {dtypes::CPn, "cp"},
    {dtypes::KPI, "kpi"},
    {dtypes::PIK, "kpi"},
    {dtypes::KsPIPI, "kspipi"},
    {dtypes::Dh, "dh"},
    {dtypes::DhCPp, "dhcp"},
    {dtypes::DhCPn, "dhcp"},
    {dtypes::DKs, "dks"},
    {dtypes::DCPpKs, "dcppks"},
    {dtypes::DCPnKs, "dcpnks"},
};

// CPV and b -> c ubar d parameters
double Cfg::rd_cp = 1.;
double Cfg::rd_kpi = 0.063;
double Cfg::deld_nega_cp = 180.;  // deg
double Cfg::deld_posi_cp = 0.;  // deg
double Cfg::deld_kpi = 10.;  // deg
double Cfg::rb = 0.02;
double Cfg::delb = 90.;  // deg
double Cfg::beta = 22.;  // deg
double Cfg::ckmgamma = 71.;  // deg
double Cfg::dtlim = 25;  // ps

double Cfg::charm_x = 0.01;
double Cfg::charm_y = 0.01;

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
//         << "mean " << mean << endl
//         << "sigma " << sigma << endl
//         << "fbkg " << fbkg << endl
//         << "wtag " << wtag << endl
         << "x " << charm_x << endl
         << "y " << charm_y << endl;
}

unique_ptr<libTatami::ToyPdf> Cfg::pdf(Cfg::ExpSetup exp) {
    const auto& expCfg = expCfgMap.at(exp);
    return std::make_unique<libTatami::ToyPdf>(expCfg.mean(), expCfg.sigma(),
                             expCfg.fbkg(), expCfg.wtag(), -dtlim, dtlim);
}

int Cfg::n_posi_cp() {return ncp / 2;}
int Cfg::n_nega_cp() {return ncp / 2;}
int Cfg::n_flv_rs() {return nfl / 2;}
int Cfg::n_flv_ws() {return nfl * rd_kpi / 2;}

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
        cout << "Cfg::nevts: Unknown eventtype "
             << static_cast<int>(type) << endl;
        return std::make_pair(0, 0);
    }
}

string Cfg::dfile(dtypes type, uint32_t nevt) {
    string output = data_path;
    switch (type) {
    case dtypes::CPp:    output += posi_cp_file; break;
    case dtypes::CPn:    output += nega_cp_file; break;
    case dtypes::KPI:    output += kpi_file; break;
    case dtypes::PIK:    output += pik_file; break;
    case dtypes::KsPIPI: output += kspp_file; break;
    default:
        cout << "Cfg::dfile: Unknown eventtype "
             << static_cast<int>(type) << endl;
        return data_path + data_file;
    }
    if (nevt)
        return output + to_string(nevt) + ".txt";
    return output + ".txt";
}

unique_ptr<DDBPars> Cfg::wfpars(const string& dcfg, const string& bcfg) {
    return make_unique<DDBPars>(rb, beta * dtr, ckmgamma * dtr,
                                delb * dtr, dcfg, bcfg);
}

unique_ptr<DDBPars> Cfg::wfpars(uint16_t seed, uint16_t idx) {
    return wfpars(get_dcfg(), get_bcfg(seed, idx));
}

unique_ptr<DDMPars> Cfg::cmpars(const string& dcfg, const string& bcfg) {
    return make_unique<DDMPars>(beta * dtr, charm_x, charm_y, dcfg, bcfg);
}

unique_ptr<DDMPars> Cfg::cmpars() {
    return cmpars(get_dcfg(), get_bcfg(8648, 0));
}

string Cfg::get_bcfg(uint16_t seed, uint16_t idx) {
    return pars_path + "wfpars_wf_tblSymABBC_bdpp_wf_seed_" +
            to_string(seed) + "_idx_" + to_string(idx) + ".txt";
}

string Cfg::get_dcfg() {
    return pars_path + dpars;
}

string Cfg::wfdtdist(uint16_t seed, uint16_t idx, double rb, dtypes type) {
    return data_path + "wf_rb_" + to_string(rb) + "_seed_"
            + to_string(seed) + "_idx_" + to_string(idx)
            + dtmap.at(type) + ".txt";
}

string Cfg::cmdtdist(double x, double y, dtypes type) {
    return data_path + "cm_x_" + to_string(x) + "_y_" + to_string(y)
            + dtmap.at(type) + ".txt";
}
