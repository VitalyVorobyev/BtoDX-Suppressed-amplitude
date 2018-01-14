#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <map>

#include "absddpars.h"

/**
 * @brief The DDMPars class manages binned parameters of double
 * Dalitz plot and accounts for the charm mixing
 */
class DDMPars : public AbsDDPars {
    // Aliases
    using ddpair = std::pair<double, double>;
    using Cache = std::map<std::string, double>;

    Cache m_cache;
    inline double cache(const std::string& name) const {
        return m_cache.at(name);
    }

    void setBeta(double x);

    bool zeroMixing() const;

    class Params {
     public:
        double Kp, Kn, C, S;
        double Kprf, Knrf, Crf, Srf;
    };
    void set_pars(int16_t bbin, int16_t dbin) const;
    mutable Params m_pars;

    ddpair coefs_dd() const;
    ddpair coefs_cp(int16_t xid) const;
    ddpair coefs_dh() const;
    ddpair coefs_dhcp(int16_t xid) const;

    ddpair calc_ud() const;
    ddpair calc_ud_dcp(int16_t xid) const;
    ddpair calc_ud_dh() const;
    ddpair calc_ud_dh_dcp(int16_t xid) const;

    double calc_f() const;
    double calc_f_dcp(int16_t xid) const;
    double calc_f_dh(int16_t xid) const;
    double calc_f_dh_dcp(int16_t xid, int16_t xih) const;

 public:
    /**
     * @brief DDMPars Constructor
     * @param beta
     * @param x
     * @param y
     * @param dcfg
     * @param bcfg
     */
    DDMPars(double beta, double x, double y,
            const std::string& dcfg, const std::string& bcfg);
    /**
     * @brief setParam
     * @param name
     * @param val
     */
    void setParam(const std::string &name, double val) override final;
    /**
     * Pair (ccoef, scoef), where ccoef (scoef) is a coefficient
     * near cos(dm*dt) (sin(dm*dt))
     */
    ddpair coefs(int16_t bbin, int16_t dbin,
                 dtypes dt=dtypes::KsPIPI) const override;
    /**
     * @brief Unnormalized coefficients U_{ij} and D_{ij}
     * @param bbin. B Dalitz plot bin number
     * @param dbin. D Dalitz plot bin number
     * @return pair of double
     */
    ddpair rawud(int16_t bbin, int16_t dbin,
                 dtypes dt=dtypes::KsPIPI) const override;
};
