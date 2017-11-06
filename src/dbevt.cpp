#include "dbevt.h"

#include <sstream>

std::string DBEvt::asStr() const {
    std::stringstream out;
    out << m_time << " " << m_tag << " " << m_bbin << " " << m_bin;
    return out.str();
}

std::ofstream& operator<< (std::ofstream& os, const DBEvt& ev) {
    os << ev.asStr();
    return os;
}

std::ostream& operator<< (std::ostream& os, const DBEvt& ev) {
    os << "Evebnt: t " << ev.m_time
       << ", tag " << ev.m_tag
       << ", bbin " << ev.m_bbin
       << ", dbin " << ev.m_bin;
    return os;
}

std::ifstream& operator>> (std::ifstream& is, DBEvt& ev) {
    is >> ev.m_time >> ev.m_tag >> ev.m_bbin >> ev.m_bin;
    return is;
}
