// test floating point exception
#pragma once
#include "../SLISC/scalar_arith.h"

void test_except()
{
    using namespace slisc;
    // test nan
    if (!ISNAN(NaN)) SLS_ERR("failed!");
    if (ISNAN(0)) SLS_ERR("failed!");
    
    Doub x;

    SLS_WARN("comment out the following lines one by one");
    x = exp(1e10); SLS_ERR("failed!"); // test overflow
    x = exp(-1e10); SLS_ERR("failed!"); // test underflow
    x = sqrt(-1.); SLS_ERR("failed!"); // test invalid
    // x = 1. / 0.; SLS_ERR("failed!"); // test division by 0
    // i = 1 / 0; SLS_ERR("failed!"); // test integer division by 0
    ++x;
}
