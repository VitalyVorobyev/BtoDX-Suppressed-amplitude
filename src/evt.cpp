#include "evt.h"

#include <sstream>

std::string Evt::asStr() const {
    std::stringstream out;
    out << m_time << " " << m_tag << " "
        << static_cast<int16_t>(m_type);
    return out.str();
}

std::ofstream& operator<< (std::ofstream& os, const Evt& ev) {
    os << ev.asStr();
    return os;
}

std::ifstream& operator>> (std::ifstream& is, Evt& ev) {
    int16_t type;
    is >> ev.m_time >> ev.m_tag >> type;
    ev.m_type = static_cast<dtypes>(type);
    return is;
}
