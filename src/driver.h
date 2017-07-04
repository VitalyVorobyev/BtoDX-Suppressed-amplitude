#ifndef DRIVER_H
#define DRIVER_H

#include <cstdint>

#include "fitmodes.h"

class Driver {
 public:
    Driver() {}
    /**
     * @brief delb_scan
     * @param type
     */
    void delb_scan(dtypes type, fitmode mode = fitmode::Simple);
    /**
     * @brief single_fit
     */
    void single_fit(dtypes type, fitmode mode = fitmode::Simple);
    /**
     * @brief cpp_and_cpn
     */
    void cpp_and_cpn(fitmode mode = fitmode::Simple);
    /**
     * @brief all_in_one. Joint fit of CP+, CP- and flavor specific events
     * @param mode
     */
    void all_in_one(fitmode mode);
    void kspp_fit(fitmode mode, uint32_t nevt, uint16_t bin=0);
    void kspp_bin_scan(fitmode mode, uint32_t nevt);
};

#endif  // DRIVER_H
