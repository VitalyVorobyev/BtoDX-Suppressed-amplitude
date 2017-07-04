#ifndef TIMEDALITZPDF_H
#define TIMEDALITZPDF_H

#include <random>       // std::default_random_engine

#include "mylibs/libTatami/toypdf.h"
#include "mylibs/libDalitz/kspipimodel.h"
#include "mylibs/libDalitz/randomdalitzpoint.h"

#include "lambda.h"

class TimeDalitzPDF {
    static constexpr double dtlo = -70;
    static constexpr double dthi =  70;

    libTatami::ToyPdf m_time_pdf;
//    KspipiModel m_dalitz_pdf;
//    Lambda_b2cud m_lambda;
//    RandomDalitzPoint m_dalitz_rndm;

    std::default_random_engine re;
    std::uniform_real_distribution<double> unif;

    auto get_dt();

public:
    TimeDalitzPDF();

    class Evt {
     public:
        double mp;
        double mm;
        double dt;

        Evt(double mp_, double mm_, double dt_) : mp(mp_), mm(mm_), dt(dt_) {}
    };
    /**
     * @brief Get PDF value for m+, m- and dt values
     */
    auto pdf(const Evt& evt);
    /**
     * @brief Generate an event accoding to the 3D PDF
     */
    auto getEvent();
};

#endif // TIMEDALITZPDF_H
