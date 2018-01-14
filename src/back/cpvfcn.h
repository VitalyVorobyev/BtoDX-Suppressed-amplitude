#ifndef CPVFCN_H
#define CPVFCN_H

#include <vector>
#include <cstdint>

#include "Minuit2/FCNBase.h"

#include "bevt.h"
#include "fitmodes.h"

class CPVFcn : public ROOT::Minuit2::FCNBase {
    const fitmode mode;
    const std::vector<BEvt>& sevtv;
    const uint16_t single_bin;

    double theErrorDef;

 public:
    CPVFcn(fitmode fmode, const std::vector<BEvt>& evtv, uint16_t sb);

    virtual double Up(void) const override {return theErrorDef;}
    virtual double operator()(const std::vector<double>& par) const override;
};

#endif // CPVFCN_H
