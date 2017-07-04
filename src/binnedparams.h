#ifndef BINNEDPARAMS_H
#define BINNEDPARAMS_H

#include <cstdint>
#include <vector>
#include <string>
#include <unordered_map>

class BinnedParams {
    const uint16_t nbins;
    double rb;
    double delb;
    double beta;
    double gamma;
    std::unordered_map<std::string, double> precl;
    std::unordered_map<int16_t, double> ccoefpcl;  // precalculated ccoefs
    std::unordered_map<int16_t, double> scoefpcl;  // precalculated scoefs
    std::unordered_map<std::string, std::vector<double>> params;

    auto scoef(int16_t bin) const;

    // internal interfaces
    inline auto get(std::string&& label, int16_t bin) const;
    inline auto prc(std::string&& label) const;

    auto ccoef(double Kp, double Kn,
               double C, double S, double C1, double S1) const;
    auto scoef(double Kp, double Kn, double Kt,
               double C, double S, double C2, double S2) const;

    bool updated;
    void update();

public:
    explicit BinnedParams(uint16_t n);

    void setBeta(double b);
    void setGamma(double b);
    void setRb(double b);
    void setDeltaB(double b);
    void setParam(std::string label, uint16_t bin, double val);
    void setParam(std::string label, std::vector<double> parv);

    std::pair<double, double> EvtFrac(int16_t bin) const;
    std::pair<double, double> coefs(int16_t bin);

    void print();
};

#endif // BINNEDPARAMS_H
