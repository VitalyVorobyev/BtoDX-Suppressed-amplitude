#ifndef WFDRIVER_H
#define WFDRIVER_H

#include <cstdint>

#include "fitmodes.h"

/** @brief High level instructions */
class WFDriver {
    WFDriver() {}

 public:
    /**
     * @brief Generates B0 -> D0 {pi+ pi-, h0}, D0 -> Ks0 pi+ pi- events
     * with charm mixing
     * @param mode - B0 decay mode
     * @param x - charm mixing parameter
     * @param y - charm mixing parameter
     * @param nevt - number of events to generate
     */
    static void gen_dd_with_charm_mix(bmode mode, double x, double y,
                                      uint64_t nevt);
    /**
     * @brief Generates B0 -> Dcp {pi+ pi-, h0} with charm mixing
     * @param mode - B0 decay mode
     * @param x - charm mixing parameter
     * @param y - charm mixing parameter
     * @param ncpp - number of events with CP = +1
     * @param ncpn - number of events with CP = -1
     */
    static void gen_cp_with_charm_mix(bmode mode, double x, double y,
                                      uint64_t ncpp, uint64_t ncpn);
    /**
     * @brief fit_with_charm_mix
     * @param mode - fit mode
     * @param x - charm mixing parameter
     * @param y - charm mixing parameter
     * @param sb - single Dalitz plot bin
     */
    static void fit_with_charm_mix(bmode mode, fitmode fmode,
                                   double x, double y, uint16_t sb=0);
    /**
     * @brief gen_dd_with_wf
     * @param nevt
     * @param idxmin
     * @param idxmax
     * @param rb
     */
    static void gen_dd_with_wf(bmode mode, double rb, uint64_t nevt,
                               uint16_t idxmin=0, uint16_t idxmax=100);
    /**
     * @brief gen_cp_with_wf
     * @param ncpp
     * @param ncnp
     * @param idxmin
     * @param idxmax
     * @param rb
     */
    static void gen_cp_with_wf(bmode mode, double rb,
                               uint64_t ncpp, uint64_t ncnp,
                               uint16_t idxmin=0, uint16_t idxmax=100);
    /**
     * @brief fit_with_wf
     * @param mode
     * @param idxmin
     * @param idxmax
     * @param rb
     * @param sb
     */
    static void fit_with_wf(bmode mode, fitmode fmode, double rb,
                            uint16_t idxmin=0, uint16_t idxmax=100,
                            uint16_t sb = 0);
};

#endif // WFDRIVER_H
