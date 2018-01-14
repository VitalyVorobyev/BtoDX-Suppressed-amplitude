#ifndef DDFCN_H
#define DDFCN_H

#include <vector>
#include <cstdint>

#include "Minuit2/FCNBase.h"

#include "dbevt.h"
#include "fitmodes.h"
#include "absddpars.h"

class DDFcn : public ROOT::Minuit2::FCNBase {
    const fitmode m_mode;
    const std::vector<DBEvt>& m_evtv;
    const uint16_t m_single_bin;
    AbsDDPars& m_pars;

    double theErrorDef;

 public:
    DDFcn(fitmode fmode, const std::vector<DBEvt>& evtv,
          uint16_t sb, AbsDDPars& pars);

    virtual double Up(void) const override {return theErrorDef;}
    virtual double operator()(const std::vector<double>& par) const override;
};

#endif // DDFCN_H
