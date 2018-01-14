#pragma once

enum class fitmode {Full, Simple, Corrected, Approx, Dh, DhCorrected, DKs};

enum class bmode {dpp, dh, DKs};

enum class dtypes : int16_t {
    CPp    =   1, // B0 -> D0 pi+ pi-, D0 -> CP+
    CPn    =  -1, // B0 -> D0 pi+ pi-, D0 -> CP-
    KPI    =   2, // B0 -> D0 pi+ pi-, D0 -> K- pi+
    PIK    =  -2, // B0 -> D0 pi+ pi-, D0 -> K+ pi-
    KsPIPI =   0, // B0 -> D0 pi+ pi-, D0 -> Ks0 pi+ pi-
    Dh     =  10, // B0 -> D0 pi0, D0 -> Ks0 pi+ pi-
    DhCPp  =  11, // B0 -> D0 pi0, D0 -> CP+
    DhCPn  = -11, // B0 -> D0 pi0, D0 -> CP-
    DKs    = 100, // B0 -> D0 Ks0, D0 -> Ks0 pi+ pi-
    DCPpKs = 101, // B0 -> D0 Ks0, D0 -> CP+
    DCPnKs =-101, // B0 -> D0 Ks0, D0 -> CP-
};
