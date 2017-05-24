#ifndef EVT_H
#define EVT_H

enum dtypes : int {CPp = +1, CPn = -1, KPI = +2, PIK = -2, KsPIPI = 0};

class Evt {
 public:
    Evt(double _time, int _tag, dtypes _type) :
        time(_time), tag(_tag), type(_type) {}
    /**
     * @brief Time difference between signal and tag B mesons
     */
    double time;
    /**
     * @brief B0 flavor tag = +1 or -1
     */
    int tag;
    /**
     * @brief D0 final state type
     * type = +1 (-1) D -> CP+ (CP-)
     * type = +2 (-2) D0 -> K-pi+ (K+pi-)
     * type = 0 -> D0 -> Ks0pi+pi-
     */
    dtypes type;
};

#endif  // EVT_H
