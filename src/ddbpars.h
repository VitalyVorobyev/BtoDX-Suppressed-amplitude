#ifndef DDBPARS_H
#define DDBPARS_H

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

#include "absddpars.h"
#include "fitmodes.h"

/**
 * @brief The DDBPars class manages binned parameters of binned double
 * Dalitz plot.
 */
class DDBPars : public AbsDDPars {
    using ddpair = std::pair<double, double>;

    double m_rb;  // rB
    double m_delb;  // average phase diff
    double m_gamma;  // CKM phase

    std::unordered_map<std::string, double> m_cache;
    void cache_upd();

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

    /** B0 -> [Ks0 pi+ pi-]_D0 pi+ pi- */
    ddpair calc_ud() const;
    double calc_f() const;
    ddpair coefs_dd() const;

    /** B0 -> Dcp pi+ pi- */
    ddpair calc_ud_cp(int xiD) const;
    double calc_f_cp(int xiD) const;
    ddpair coefs_cp(int xiD) const;

    /** B0 -> D0 h0, D0 -> Ks0 pi+ pi- */
    ddpair calc_ud_dh() const;
    double calc_f_dh() const;
    ddpair coefs_dh() const;

    /** B0 -> Dcp h0 */
    ddpair calc_ud_dhcp(int xiD) const;
    double calc_f_dhcp(int xiD) const;
    ddpair coefs_dhcp(int xiD) const;

    mutable Params m_pars;

public:
    DDBPars(double rb, double beta, double gamma, double delb,
            const std::string& dcfg, const std::string& bcfg);

    void set_beta(double x) override;
    void set_rb(double x);
    void set_delb(double x);
    void set_gamma(double x);
    void set(double rb, double _beta, double gamma);

    double rb() const {return m_rb;}
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

#endif  // DDBPARS_H
