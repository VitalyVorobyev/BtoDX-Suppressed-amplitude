#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <map>

#include "absddpars.h"
#include "fitmodes.h"

/** @brief The DDBPars class manages binned parameters of
 *  binned double Dalitz plot. */
class DDBPars : public AbsDDPars {
    using ddpair = std::pair<double, double>;
    using Cache = std::map<std::string, double>;

    Cache m_cache;
    void updateCache();
    inline double cache(const std::string& name) const {
        return m_cache.at(name);
    }

    class Params {
     public:
        Params();
        double Kp, Kn, C, S;
        double Kprf, Knrf, Kpwf, Knwf;
        double Crf, Srf, Cwf, Swf;
        double Ctp, Stp, Ctn, Stn;
        double Cpp, Spp, Cpn, Spn;
    };

    void set_pars(int16_t bbin, int16_t dbin) const;

    /** Auxiliary methods */
    void setRb(double x);
    void setGamma(double x);
    void setBeta(double x);
    void setDeltaB(double x);

    /** B0 -> [D0 -> Ks0 pi+ pi-] pi+ pi- */
    ddpair calc_ud() const;
    double calc_f() const;
    ddpair coefs_dd() const;

    /** B0 -> Dcp pi+ pi- */
    ddpair calc_ud_cp(int xiD) const;
    double calc_f_cp(int xiD) const;
    ddpair coefs_cp(int xiD) const;

    /** B0 -> [D0 -> Ks0 pi+ pi-] {h0, Ks0} */
    ddpair calc_ud_dh() const;
    double calc_f_dh() const;
    ddpair coefs_dh() const;

    /** B0 -> Dcp {h0, Ks0} */
    ddpair calc_ud_dhcp(int xiD) const;
    double calc_f_dhcp(int xiD) const;
    ddpair coefs_dhcp(int xiD) const;

    mutable Params m_pars;

 public:
    /**
     * @brief DDBPars Constructor
     * @param rb
     * @param beta
     * @param gamma
     * @param delb
     * @param dcfg
     * @param bcfg
     */
    DDBPars(double rb, double beta, double gamma, double delb,
            const std::string& dcfg, const std::string& bcfg);

    void setParam(const std::string& name, double val) override final;
    /**
     * Pair (ccoef, scoef), where ccoef (scoef) is a coefficient
     * near cos(dm*dt) (sin(dm*dt))
     */
    ddpair coefs(int16_t bbin, int16_t dbin,
                 dtypes dt=dtypes::KsPIPI) const override final;
    /**
     * @brief Unnormalized coefficients U_{ij} and D_{ij}
     * @param bbin. B Dalitz plot bin number
     * @param dbin. D Dalitz plot bin number
     * @return pair of double
     */
    ddpair rawud(int16_t bbin, int16_t dbin,
                 dtypes dt=dtypes::KsPIPI) const override final;
};
