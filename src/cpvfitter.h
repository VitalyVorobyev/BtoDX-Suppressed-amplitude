#ifndef CPVFITTER_H
#define CPVFITTER_H

#include <vector>
#include <string>
#include <cstdint>

#include "Minuit2/MnUserParameterState.h"

#include "dbevt.h"
#include "fitmodes.h"
#include "absddpars.h"

class CPVFitter {
    std::vector<BEvt> sevtv;
    fitmode sfmode;
    uint16_t single_bin;

 public:
    explicit CPVFitter(fitmode mode = fitmode::Simple);
    CPVFitter(std::string& input_file_name, fitmode mode = fitmode::Simple);
    CPVFitter(std::vector<std::string>& input_files,
              fitmode mode = fitmode::Simple);
    /**
     * @brief ReadData. Read event vector from text file
     * @param input_file_name. Input text file
     */
    void ReadBEvt(std::string& input_file_name);
    ROOT::Minuit2::MnUserParameterState MakeFit();

    auto ReadDBEvt(const std::string& fname);
    ROOT::Minuit2::MnUserParameterState DDFit(const std::string& fname,
                                              AbsDDPars& pars);

    void setSingleBin(uint16_t sb) {single_bin = sb;}
    void unsetSingleBin() {setSingleBin(0);}

    std::pair<std::vector<double>, std::vector<double>> Scan(
            double lo=0., double hi=180., uint16_t nbin=180) const;
};

#endif  // CPVFITTER_H
