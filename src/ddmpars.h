#ifndef DDMPARS_H
#define DDMPARS_H

#include <cstdint>
#include <vector>
#include <string>

#include "absddpars.h"

/**
 * @brief The DDMPars class manages binned parameters of double
 * Dalitz plot and accounts for the charm mixing
 */
class DDMPars : public AbsDDPars {
    double m_x;  // charm mixing
    double m_y;  // charm mixing
    double m_cosdbeta;
    double m_sindbeta;
    bool m_zeromix;

    class Params {
     public:
        double Kp, Kn, C, S;
        double Kprf, Knrf, Crf, Srf;
    };
    void set_pars(int16_t bbin, int16_t dbin) const;
    mutable Params m_pars;

    std::pair<double, double> coefs_dd() const;
    std::pair<double, double> coefs_cp(int16_t xid) const;
    std::pair<double, double> coefs_dh() const;
    std::pair<double, double> coefs_dhcp(int16_t xid) const;

    std::pair<double, double> calc_ud() const;
    std::pair<double, double> calc_ud_dcp(int16_t xid) const;
    std::pair<double, double> calc_ud_dh() const;
    std::pair<double, double> calc_ud_dh_dcp(int16_t xid) const;

    double calc_f() const;
    double calc_f_dcp(int16_t xid) const;
    double calc_f_dh(int16_t xid) const;
    double calc_f_dh_dcp(int16_t xid, int16_t xih) const;

 public:
    DDMPars(double beta, double x, double y,
            const std::string& dcfg, const std::string& bcfg);
    void set_xy(double x, double y);
    void set_beta(double x);
    /**
     * Pair (ccoef, scoef), where ccoef (scoef) is a coefficient
     * near cos(dm*dt) (sin(dm*dt))
     */
    std::pair<double, double> coefs(int16_t bbin, int16_t dbin,
                                    dtypes dt=dtypes::KsPIPI) const override;
    /**
     * @brief Unnormalized coefficients U_{ij} and D_{ij}
     * @param bbin. B Dalitz plot bin number
     * @param dbin. D Dalitz plot bin number
     * @return pair of double
     */
    std::pair<double, double> rawud(int16_t bbin, int16_t dbin,
                                    dtypes dt=dtypes::KsPIPI) const override;
};

#endif // DDMPARS_H
