#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <map>

#include "fitmodes.h"

/**
 * @brief The DDMPars class manages binned parameters of double Dalitz plot
 */
class AbsDDPars {
    // Aliases
    using ParamsMap = std::map<std::string, double>;
    using ParsVecMap = std::map<std::string, std::vector<double>>;
    /** Number of bins on B meson DP */
    static constexpr uint16_t m_nbbins = 8;
    /** Number of bins on D meson DP */
    static constexpr uint16_t m_ndbins = 8;
    /** CKM phase */
    ParamsMap m_pars;
    /** Bin integrals vectors */
    ParsVecMap m_int;

    int16_t readDConfig(const std::string& fname, bool verb);
    int16_t readBConfig(const std::string& fname, bool verb);

 public:
    AbsDDPars(const std::string& dcfg, const std::string& bcfg,
              double beta=22.);
    virtual ~AbsDDPars() = default;

    void set_c(std::vector<double>& x);
    void set_s(std::vector<double>& x);
    const std::vector<double>& get_c() const;
    const std::vector<double>& get_s() const;

    virtual void setParam(const std::string& name, double val);
    void setNewParam(const std::string& name, double val);
    double getParam(const std::string& name) const;
    double getInt(const std::string& name, int16_t bin) const;
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
