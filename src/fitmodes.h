#ifndef FITMODES_H
#define FITMODES_H

enum class fitmode {Full, Simple, Corrected, Approx, Dh, DhCorrected};

enum class bmode {dpp, dh};

enum class dtypes : int16_t {
    CPp    =  1, // B0 -> D0 pi+ pi-, D0 -> CP+
    CPn    = -1, // B0 -> D0 pi+ pi-, D0 -> CP-
    KPI    =  2, // B0 -> D0 pi+ pi-, D0 -> K- pi+
    PIK    = -2, // B0 -> D0 pi+ pi-, D0 -> K+ pi-
    KsPIPI =  0, // B0 -> D0 pi+ pi-, D0 -> Ks0 pi+ pi-
    Dh     = 10, // B0 -> D0 pi0, D0 -> Ks0 pi+ pi-
    DhCPp  = 11, // B0 -> D0 pi0, D0 -> CP+
    DhCPn  =-11  // B0 -> D0 pi0, D0 -> CP-
};

#endif // FITMODES_H
