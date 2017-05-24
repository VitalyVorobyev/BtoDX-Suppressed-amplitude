#ifndef CPVGEN_H
#define CPVGEN_H

#include <utility>  // pair
#include <vector>
#include <string>

#include "mylibs/libTatami/toypdf.h"

#include "cfg.h"
#include "evt.h"

/**
 * @brief The CPVGen class.
 */
class CPVGen {
 public:
    CPVGen() {}
    /**
     * @brief Generate vector of Evt events with
     * custum configuration
     * @param pdf
     * @param nev
     * @param evtv
     * @param type. Type of D meson final state
     */
    void Gen(libTatami::ToyPdf& pdf, std::pair<int, int> nev,
             std::vector<Evt>& evtv, dtypes type);
    /**
     * @brief Generate vector of Evt events with
     * predefined configuration
     * @param evtv
     * @param type. Type of D meson final state
     */
    void Gen(std::vector<Evt>& evtv, dtypes type);
    /**
     * @brief Write events into text file in random order
     * @param evts. Vector of events
     * @param fname. File name
     */
    void ShuffleAndWrite(std::vector<Evt>& evts, std::string& fname);
    /**
     * @brief Generate events with custum configuration and
     * save it in file
     * @param pdf PDF
     * @param nev - pair of integers specifying number of B0 and B0-bar
     * tagged events
     * @param fname - output file name
     * @param type - type of D meson final state
     */
    void GenAndWrite(libTatami::ToyPdf& pdf, std::pair<int, int> nev,
                     std::string& fname, dtypes type);
    /**
     * @brief Generate events with predefined configuration and
     * save it in file
     * @param type - type of D meson final state
     */
    void GenAndWrite(dtypes type);
    /**
     * @brief GenAndWrite
     * @param type - type of D meson final state
     * @param nev - number of events
     */
    void GenAndWrite(dtypes type, std::pair<int, int> nev);
    /**
     * @brief Generate all types of events and save into files
     */
    void CompleteData();
};

#endif  // CPVGEN_H
