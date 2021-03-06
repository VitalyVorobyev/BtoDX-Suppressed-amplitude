#include "wfdriver.h"

#include <iostream>
#include <map>
#include <memory>

//#include "cpvgen.h"
#include "ddgev.h"
#include "cfg.h"
//#include "wfminuitfitter.h"
#include "cpvminimizer.h"

#include "ddbpars.h"
#include "ddmpars.h"

using std::cout;
using std::endl;
using std::map;
using std::unique_ptr;

const map<bmode, dtypes> dtkspp = {
    {bmode::dh, dtypes::Dh},
    {bmode::dpp, dtypes::KsPIPI},
    {bmode::DKs, dtypes::DKs}
};

const map<bmode, dtypes> dtcp = {
    {bmode::dh,  dtypes::DhCPp},
    {bmode::dpp, dtypes::CPp},
    {bmode::DKs, dtypes::DCPpKs}
};

const map<fitmode, bool> zerofitrb = {
    {fitmode::Simple,      true},
    {fitmode::Full,        true},
    {fitmode::Dh,          true},
    {fitmode::Corrected,   false},
    {fitmode::DhCorrected, false},
    {fitmode::DKs,         false}
};

void WFDriver::gen_dd_with_charm_mix(bmode mode, double x, double y,
                                     uint64_t nevt) {
    Cfg::set_charm_mix(x, y);
    auto pars = Cfg::cmpars();
    auto ofile = Cfg::cmdtdist(x, y, dtkspp.at(mode));
    if (mode == bmode::dpp)
        DDGev::DoubleDalitz(nevt, *pars, ofile);
    else
        DDGev::DhDalitz(nevt, *pars, ofile);
}

void WFDriver::gen_cp_with_charm_mix(bmode mode, double x, double y,
                                     uint64_t ncpp, uint64_t ncpn) {
    Cfg::set_charm_mix(x, y);
    auto pars = Cfg::cmpars();
    auto ofile = Cfg::cmdtdist(x, y, dtcp.at(mode));
    if (mode == bmode::dpp)
        DDGev::CPDalitz(ncpp, ncpn, *pars, ofile);
    else
        DDGev::DhCP(ncpp, ncpn, *pars, ofile);
}

void WFDriver::gen_cp_with_wf(bmode mode, double rb,
                              uint64_t ncpp, uint64_t ncpn,
                              uint16_t idxmin, uint16_t idxmax) {
    cout << "gen_cp_with_wf: using rB = " << rb << endl;
    Cfg::set_rb(rb);
    if (mode == bmode::dh || mode == bmode::DKs) {
        idxmin = 0;
        idxmax = 1;
    }
    for (auto idx = idxmin; idx < idxmax; idx++) {
        cout << "idx " << idx + 1 << " / " << idxmax << endl;
        auto pars = Cfg::wfpars(8648, idx);
        auto ofile = Cfg::wfdtdist(8648, idx, rb, dtcp.at(mode));
        cout << "Events to be written in " << ofile << endl;
        if (mode == bmode::dpp)
            DDGev::CPDalitz(ncpp, ncpn, *pars, ofile);
        else
            DDGev::DhCP(ncpp, ncpn, *pars, ofile);
    }
}

void WFDriver::gen_dd_with_wf(bmode mode, double rb, uint64_t nevt,
                              uint16_t idxmin, uint16_t idxmax) {
    cout << "gen_dd_with_wf: using rB = " << rb << endl;
    Cfg::set_rb(rb);
    if (mode == bmode::dh || mode == bmode::DKs) {
        idxmin = 0;
        idxmax = 1;
    }
    for (auto idx = idxmin; idx < idxmax; idx++) {
        cout << "idx " << idx + 1 << " / " << idxmax << endl;
        auto pars = Cfg::wfpars(8648, idx);
        auto ofile = Cfg::wfdtdist(8648, idx, rb, dtkspp.at(mode));
        cout << ofile << endl;
        if (mode == bmode::dpp)
            DDGev::DoubleDalitz(nevt, *pars, ofile);
        else
            DDGev::DhDalitz(nevt, *pars, ofile);
    }
}

void WFDriver::fit_with_charm_mix(bmode mode, fitmode fmode,
                                  double x, double y, uint16_t sb) {
    if (fmode == fitmode::Simple || fmode == fitmode::Full || fmode == fitmode::Dh) {
        cout << "fit_with_charm_mix: zero mix" << endl;
        Cfg::set_charm_mix(0., 0.);
    } else if (fmode == fitmode::Corrected || fmode == fitmode::DhCorrected) {
        cout << "fit_with_charm_mix: Corrected" << endl;
        Cfg::set_charm_mix(x, y);
    }
    using Fitter = CPVMinimizer;
    Fitter::FlushData();
    unique_ptr<AbsDDPars> pars = Cfg::cmpars();
    Fitter::ReadData(Cfg::cmdtdist(x, y, dtcp.at(mode)));
    if (mode == bmode::dpp)
        Fitter::ReadBData(Cfg::cmdtdist(x, y, dtkspp.at(mode)));
    else
        Fitter::ReadData(Cfg::cmdtdist(x, y, dtkspp.at(mode)));
    Fitter::MakeFit(pars, fmode, sb);
}

void WFDriver::fit_with_wf(bmode mode, fitmode fmode, double rb,
                           uint16_t idxmin, uint16_t idxmax, uint16_t sb) {
    if (zerofitrb.at(fmode)) {
        cout << "fit_with_wf: zero rb" << endl;
        Cfg::set_rb(0.);
    } else {
        cout << "fit_with_wf: using rb = " << rb << endl;
        Cfg::set_rb(rb);
    }
    if (mode == bmode::dh || mode == bmode::DKs) {
        idxmin = 0;
        idxmax = 1;
    }
    using Fitter = CPVMinimizer;
    for (auto idx = idxmin; idx < idxmax; idx++) {
        cout << "idx " << idx + 1 << " / " << idxmax << endl;
        Fitter::FlushData();
        unique_ptr<AbsDDPars> pars = Cfg::wfpars(8648, idx);
        Fitter::ReadData(Cfg::wfdtdist(8648, idx, rb, dtcp.at(mode)));
        if (mode == bmode::dpp)
            Fitter::ReadBData(Cfg::wfdtdist(8648, idx, rb, dtkspp.at(mode)));
        else
            Fitter::ReadData(Cfg::wfdtdist(8648, idx, rb, dtkspp.at(mode)));
        Fitter::MakeFit(pars, fmode, sb);
    }
}
