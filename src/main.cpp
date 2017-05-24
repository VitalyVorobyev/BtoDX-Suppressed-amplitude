#include "driver.h"

int main(int argc, char** argv) {
    Cfg::print_config();
    Driver driver;
//    driver.delb_scan(dtypes::CPp);
//    driver.single_fit(dtypes::CPp);
//    driver.single_fit(dtypes::CPn);
//    driver.single_fit(dtypes::KPI, fitmode::Corrected);
//    driver.single_fit(dtypes::PIK);
//    driver.cpp_and_cpn(fitmode::Simple);
//    driver.cpp_and_cpn(fitmode::Corrected);
//    driver.cpp_and_cpn(fitmode::Full);
//    driver.all_in_one(fitmode::Full);
    driver.all_in_one(fitmode::Simple);
    return 0;
}
