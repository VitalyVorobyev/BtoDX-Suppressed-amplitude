#include "bevt.h"

#include <sstream>

std::string BEvt::asStr() const {
    std::stringstream out;
    out << m_time << " " << m_tag << " " << m_bin << " "
        << static_cast<int16_t>(m_type);
    return out.str();
}

std::ofstream& operator<< (std::ofstream& os, const BEvt& ev) {
    os << ev.asStr();
    return os;
}

std::ostream& operator<< (std::ostream& os, const BEvt& ev) {
    os << "Event: t " << ev.m_time
       << ", tag " << ev.m_tag
       << ", bin " << ev.m_bin
       << ", type " << static_cast<int>(ev.m_type);
    return os;
}

std::ifstream& operator>> (std::ifstream& is, BEvt& ev) {
    int16_t type;
    is >> ev.m_time >> ev.m_tag >> ev.m_bin >> type;
    ev.m_type = static_cast<dtypes>(type);
    return is;
}
