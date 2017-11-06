#ifndef EVT_H
#define EVT_H

#include <iostream>
#include <cstdint>
#include <string>
#include <memory>
#include <vector>
#include <fstream>

#include "fitmodes.h"

/**
 * @brief The event class for time-dependent CPV fit
 */
class Evt {
 protected:
    /**
     * @brief Time difference between signal and tag B mesons
     */
    double m_time;
    /**
     * @brief B0 flavor tag = +1 or -1
     */
    int16_t m_tag;
    /**
     * @brief D0 final state type
     * type = +1 (-1) : D0 -> CP+ (CP-)
     * type = +2 (-2) : D0 -> K-pi+ (K+pi-)
     * type =  0      : D0 -> Ks0pi+pi-
     */
    dtypes m_type;

 public:
    /** Constructor */
    Evt(double time, int16_t tag, dtypes type) :
        m_time(time), m_tag(tag), m_type(type) {}
    /** @brief Time */
    auto t() const {return m_time;}
    /** @brief initial B meson flavor */
    auto Tag() const {return m_tag;}
    /** @brief D0 decay mode */
    auto Type() const {return m_type;}
    /** Event as text string */
    virtual std::string asStr() const;

    friend std::ofstream& operator <<(std::ofstream& os, const Evt& ev);
    friend std::ifstream& operator >>(std::ifstream& is, Evt& ev);
};

template<class T>
std::vector<T> ReadEvents(const std::string& ifname) {
    std::vector<T> events;
    std::ifstream ifile(ifname, std::ifstream::in);
    if (!ifile.is_open()) {
        std::cerr << "Can't open file " << ifname << std::endl;
        return events;
    }
    std::cout << "ReadEvents: read data from " << ifname << std::endl;
    T evt;
    while (ifile.good()) {
        ifile >> evt;
        events.emplace_back(evt);
    }
    std::cout << events.size() << " events found" << std::endl;
    ifile.close();
    return events;
}

template<class T>
std::vector<T> ReadEvents(const std::vector<std::string>& ifnames) {
    std::vector<T> events;
    for (const auto& fname : ifnames) {
        auto evv = ReadEvents<T>(fname);
        events.insert(events.end(), evv.beind(), evv.end());
    }
    return events;
}

#endif  // EVT_H
