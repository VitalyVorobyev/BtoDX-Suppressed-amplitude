#ifndef CPVGEN_H
#define CPVGEN_H

#include <utility>  // pair
#include <vector>
#include <string>

#include "mylibs/libTatami/toypdf.h"

#include "bevt.h"

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
    void Gen(libTatami::ToyPdf& pdf, std::pair<uint64_t, uint64_t> &&nev,
             std::vector<BEvt>& evtv, dtypes type, int16_t bin=0);
    /**
     * @brief Generate vector of Evt events with
     * predefined configuration
     * @param evtv
     * @param type. Type of D meson final state
     */
    void Gen(std::vector<BEvt>& evtv, dtypes type);
    /**
     * @brief Write events into text file in random order
     * @param evts. Vector of events
     * @param fname. File name
     */
    void ShuffleAndWrite(std::vector<BEvt>& evts, const std::string& fname);
    /**
     * @brief Generate events with custum configuration and
     * save it in file
     * @param pdf PDF
     * @param nev - pair of integers specifying number of B0 and B0-bar
     * tagged events
     * @param fname - output file name
     * @param type - type of D meson final state
     */
    void GenAndWrite(libTatami::ToyPdf& pdf,
                     std::pair<uint64_t, uint64_t>&& nev,
                     const std::string &fname, dtypes type);
    /**
     * @brief Generate events with predefined configuration and
     * save it in file
     * @param type - type of D meson final state
     */
    void GenAndWrite(dtypes type, const std::string &pref=default_str);
    /**
     * @brief GenAndWrite
     * @param type - type of D meson final state
     * @param nev - number of events
     */
    void GenAndWrite(dtypes type, std::pair<uint64_t, uint64_t>&& nev,
                     const std::string& fname=default_str);
    /**
     * @brief GenKspipi
     * @param nev
     * @param fname
     */
    void GenKspipi(std::pair<uint64_t, uint64_t>&& nev,
                   const std::string& fname=default_str);
    /**
     * @brief Generate all types of events and save into files
     */
    void CompleteData(const std::string &pref=default_str);

 private:
    static const std::string default_str;
};

#endif  // CPVGEN_H
