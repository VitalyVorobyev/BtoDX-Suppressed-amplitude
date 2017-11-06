#ifndef DBEVT_H
#define DBEVT_H

#include "bevt.h"
#include <fstream>
#include <iostream>

/**
 * @brief The DBEvt class. Toy MC event for B0 -> D0 pi+ pi-, D0 -> Ks0 pi+ pi-
 */
class DBEvt : public BEvt {
    /** @brief B0 meson Dalitz plot bin number */
    int16_t m_bbin;

public:
    DBEvt(double time=0, int16_t tag=0, int16_t dbin=0, int16_t bbin=0) :
        BEvt(time, tag, dbin, dtypes::KsPIPI), m_bbin(bbin) {}
    /** @brief B0 meson Dalitz plot bin number */
    auto Bbin() const {return m_bbin;}
    /** @brief D0 meson Dalitz plot bin number */
    auto Dbin() const {return Bin();}
    /** Event as text string */
    std::string asStr() const override;

    friend std::ofstream& operator<< (std::ofstream& os, const DBEvt& ev);
    friend std::ostream& operator<< (std::ostream& os, const DBEvt& ev);
    friend std::ifstream& operator>> (std::ifstream& is, DBEvt& ev);
};

#endif // DBEVT_H
