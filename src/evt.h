#ifndef EVT_H
#define EVT_H

#include <cstdint>

#include "fitmodes.h"

/**
 * @brief The event class for time-dependent CPV fit
 */
class Evt {
 public:
    Evt(double _time, int16_t _tag, dtypes _type) :
        time(_time), tag(_tag), type(_type) {}
    /**
     * @brief Time difference between signal and tag B mesons
     */
    double time;
    /**
     * @brief B0 flavor tag = +1 or -1
     */
    int16_t tag;
    /**
     * @brief D0 final state type
     * type = +1 (-1) D -> CP+ (CP-)
     * type = +2 (-2) D0 -> K-pi+ (K+pi-)
     * type = 0 -> D0 -> Ks0pi+pi-
     */
    dtypes type;
};

#endif  // EVT_H
