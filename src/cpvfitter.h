#ifndef CPVFITTER_H
#define CPVFITTER_H

#include <vector>
#include <string>
#include <cstdint>

#include "Minuit2/MnUserParameterState.h"

#include "bevt.h"
#include "fitmodes.h"

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
    void ReadData(std::string& input_file_name);
    /**
     * @brief TheFCN. Unbinned Log Likelihood Function
     * @param npar number of parameters to fit
     * @param grad precomputed gradient (optional)
     * @param fval fcn value (to be filled)
     * @param p parameters vector
     * @param iflag flag
     */
    ROOT::Minuit2::MnUserParameterState MakeFit();

    void setSingleBin(uint16_t sb) {single_bin = sb;}
    void unsetSingleBin() {setSingleBin(0);}
};

#endif  // CPVFITTER_H
