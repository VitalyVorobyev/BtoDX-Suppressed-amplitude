#ifndef ABSDDPARS_H
#define ABSDDPARS_H

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "fitmodes.h"

/**
 * @brief The DDMPars class manages binned parameters of double Dalitz plot
 */
class AbsDDPars {
    /** Number of bins on B meson DP */
    static constexpr uint16_t m_nbbins = 8;
    /** Number of bins on D meson DP */
    static constexpr uint16_t m_ndbins = 8;
    /** CKM phase */
    double m_beta;

    int16_t readDConfig(const std::string& fname, bool verb);
    int16_t readBConfig(const std::string& fname, bool verb);

 protected:
    /** Bin integrals vectors */
    std::unordered_map<std::string, std::vector<double>> m_int;

 public:
    AbsDDPars(const std::string& dcfg, const std::string& bcfg);

    virtual void set_beta(double x) {m_beta = x;}
    double beta() const {return m_beta;}
    void set_c(std::vector<double>& x);
    void set_s(std::vector<double>& x);
    const std::vector<double>& get_c() const;
    const std::vector<double>& get_s() const;
    /**
     * Pair (ccoef, scoef), where ccoef (scoef) is a coefficient
     * near cos(dm*dt) (sin(dm*dt))
     */
    virtual std::pair<double, double> coefs(int16_t bbin, int16_t dbin,
                                            dtypes dt=dtypes::KsPIPI) const = 0;
    /**
     * @brief Unnormalized coefficients U_{ij} and D_{ij}
     * @param bbin. B Dalitz plot bin number
     * @param dbin. D Dalitz plot bin number
     * @return pair of double
     */
    virtual std::pair<double, double> rawud(int16_t bbin, int16_t dbin,
                                            dtypes dt=dtypes::KsPIPI) const = 0;
    void print() const;

    static uint16_t nbbins() {return m_nbbins;}
    static uint16_t ndbins() {return m_ndbins;}
};

#endif // ABSDDPARS_H
