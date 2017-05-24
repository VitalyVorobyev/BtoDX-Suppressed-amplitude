#ifndef CPVFITTER_H
#define CPVFITTER_H

#include <vector>
#include <string>
#include <unordered_map>

#include "mylibs/libTatami/toypdf.h"

#include "lambda.h"
#include "evt.h"
#include "cfg.h"

enum fitmode : int {Full, Simple, Corrected};

class CPVFitter {
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
    void MakeFit();

 private:
    fitmode fmode;
    libTatami::ToyPdf pdf;
    std::vector<Evt> evtv;
    std::unordered_map<dtypes, Lambda_b2cud> lambdas;

    // Have to make static interface for TVirtualFitter
    static void TheFCN(int& npar, double* grad, double &fval,
                       double *p, int iflag);
    static libTatami::ToyPdf* spdf;
    static std::vector<Evt>* sevtv;
    static std::unordered_map<dtypes, Lambda_b2cud>* slambdas;
    static fitmode sfmode;
};

#endif  // CPVFITTER_H
