//#include "driver.h"

//#include <iostream>
//#include <string>
//#include <utility>  // pair
//#include <cmath>

//#include "cpvgen.h"
//#include "ddgev.h"
//#include "cpvfitter.h"
//#include "evt.h"
//#include "cfg.h"
//#include "minuitfitter.h"

//#include "ddbpars.h"
//#include "ddmpars.h"

//using std::string;
//using std::vector;

//using std::cout;
//using std::endl;

//using std::make_pair;

//string mode_lbl(fitmode mode) {
//    if (mode == fitmode::Full) return "1";
//    if (mode == fitmode::Simple) return "2";
//    if (mode == fitmode::Corrected) return "3";
//    return "unknown_mode";
//}

//void Driver::delb_scan(dtypes type, fitmode mode) {
//    cout << "delb_scan" << endl;
//    uint64_t N = pow(10, 6);
//    auto fname = Cfg::dfile(type) + "scan" + mode_lbl(mode);
//    CPVGen gen;
//    for (double db = 0.; db < 360.; db += 15.) {
//        std::cout << "db " << db << std::endl;
//        Cfg::set_delb(db);
//        gen.GenAndWrite(type, make_pair(N / 2, N / 2), fname);
//        CPVFitter fitter(fname, mode);
//    }
//}

//void Driver::single_fit(dtypes type, fitmode mode) {
//    auto lambda = Cfg::lamf(type);
//    auto pdf = Cfg::pdf();
//    pdf.SetS(lambda.scoef());
//    pdf.SetC(lambda.ccoef());
//    const int N = pow(10, 6);
//    CPVGen gen;
//    auto fname = Cfg::dfile(type);
//    auto nevts = make_pair(N / 2, N / 2);
//    gen.GenAndWrite(pdf, nevts, fname, type);
//    CPVFitter fitter(fname, mode);
//    fitter.MakeFit();
//}

//void Driver::cpp_and_cpn(fitmode mode) {
//    CPVGen gen;
//    const int N = pow(10, 6);
//    auto nevts = make_pair(N / 2, N / 2);
//    gen.GenAndWrite(dtypes::CPp, nevts);
//    gen.GenAndWrite(dtypes::CPn, nevts);
//    std::vector<std::string> files = {
//        Cfg::dfile(dtypes::CPp),
//        Cfg::dfile(dtypes::CPn)
//    };
//    CPVFitter fitter(files, mode);
//    fitter.MakeFit();
//}

//void Driver::all_in_one(fitmode mode) {
////    Cfg::set_ncp(2*std::pow(10, 6));
////    Cfg::set_nfl(2*std::pow(10, 6));
//    CPVGen gen;
//    auto pref = mode_lbl(mode);
//    cout << pref << endl;
//    gen.CompleteData(pref);
//    vector<string> files = {
//        Cfg::dfile(dtypes::CPp) + pref,
//        Cfg::dfile(dtypes::CPn) + pref,
//        Cfg::dfile(dtypes::KPI) + pref,
//        Cfg::dfile(dtypes::PIK) + pref
//    };
//    CPVFitter fitter(files, mode);
//    fitter.MakeFit();
//}

//void Driver::kspp_fit(
//        fitmode mode, uint32_t nevt, uint16_t bin) {
//    auto fname = Cfg::dfile(dtypes::KsPIPI, nevt);
//    MinuitFitter<DDBPars> fitter(fname, mode);
//    fitter.setSingleBin(bin);
//    fitter.MakeFit();
//}

//void Driver::kspp_bin_scan(fitmode mode, uint32_t nevt) {
//    auto fname = Cfg::dfile(dtypes::KsPIPI, nevt);
//    cout << "Fit scan start" << endl;
//    for (auto bin = 1; bin <= 8; bin++) {
//        cout << "### bin " << bin << " ###" << endl;
//        CPVFitter fitter(fname, mode);
//        fitter.setSingleBin(bin);
//        fitter.MakeFit();
//    }
//}

//void Driver::kspp_bin_profile(fitmode mode, uint32_t nevt) {
//    auto fname = Cfg::dfile(dtypes::KsPIPI, nevt);
//    CPVFitter fitter(fname, mode);
//    for (auto bin = 1; bin <= 8; bin++) {
//        cout << "### Profile bin " << bin << " ###" << endl;
//        fitter.setSingleBin(bin);
//        const auto scan = fitter.Scan();
//        for (auto idx = 0; idx <= 180; idx++) cout << scan.second[idx] << ", ";
//        cout << endl;
//    }
//}

//void Driver::kspp_profile(fitmode mode, uint32_t nevt, uint16_t bin) {
//    cout << "### Profile ###" << endl;
//    auto fname = Cfg::dfile(dtypes::KsPIPI, nevt);
//    CPVFitter fitter(fname, mode);
//    fitter.setSingleBin(bin);
//    const auto scan = fitter.Scan();
//    for (auto idx = 0; idx <= 180; idx++) cout << scan.second[idx] << ", ";
//    cout << endl;
//}
