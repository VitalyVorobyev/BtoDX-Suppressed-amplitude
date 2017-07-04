#ifndef MINUITFITTER_H
#define MINUITFITTER_H

#include <vector>
#include <string>
#include <cstdint>

#include "bevt.h"
#include "fitmodes.h"

class MinuitFitter {
    static std::vector<BEvt> sevtv;
    static uint16_t single_bin;
    static fitmode sfmode;
    static void TheFCN(int& npar, double* grad,
                       double &fval, double *p, int iflag);

public:
    explicit MinuitFitter(fitmode mode = fitmode::Simple);
    MinuitFitter(const std::string& input_file_name, fitmode mode = fitmode::Simple);
    MinuitFitter(const std::vector<std::string>& input_files,
              fitmode mode = fitmode::Simple);
    /**
     * @brief ReadData. Read event vector from text file
     * @param input_file_name. Input text file
     */
    void ReadData(const std::string &input_file_name);
    /**
     * @brief TheFCN. Unbinned Log Likelihood Function
     * @param npar number of parameters to fit
     * @param grad precomputed gradient (optional)
     * @param fval fcn value (to be filled)
     * @param p parameters vector
     * @param iflag flag
     */
    void MakeFit();

    void setSingleBin(uint16_t sb) {single_bin = sb;}
    void unsetSingleBin() {setSingleBin(0);}
};

#endif // MINUITFITTER_H
