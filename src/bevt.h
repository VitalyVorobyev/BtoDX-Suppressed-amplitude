#ifndef BEVT_H
#define BEVT_H

#include "evt.h"
#include <fstream>
#include <iostream>

/**
 * @brief The event class for binned time-dependent CPV fit
 */
class BEvt : public Evt {
 protected:
    /**
     * @brief bin. If dtypes == KsPIPI, then bin is D0 Dalitz plot bin number
     * else bin is B0 Dalitz plot bin number
     */
    int16_t m_bin;

 public:
    BEvt(double time=0., int16_t tag=0, int16_t bin=0,
         dtypes type=dtypes::KsPIPI) :
        Evt(time, tag, type), m_bin(bin) {}
    /** @brief Dalitz plot bin number */
    auto Bin() const {return m_bin;}
    /** Event as text string */
    virtual std::string asStr() const override;
    /** Write to a text file */
    friend std::ofstream& operator<< (std::ofstream& os, const BEvt& ev);
    /** Print to stdout */
    friend std::ostream& operator<< (std::ostream& os, const BEvt& ev);
    /** Read from a text file */
    friend std::ifstream& operator>> (std::ifstream& is, BEvt& ev);
};

#endif // BEVT_H
