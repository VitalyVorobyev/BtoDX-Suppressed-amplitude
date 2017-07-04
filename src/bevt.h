#ifndef BEVT_H
#define BEVT_H

#include "evt.h"

/**
 * @brief The event class for binned time-dependent
 * CPV fit
 */
class BEvt : public Evt {
public:
    BEvt(double _time, int16_t _tag, int16_t _bin, dtypes _type) :
        Evt(_time, _tag, _type), bin(_bin) {}

    int16_t bin;
};

#endif // BEVT_H
