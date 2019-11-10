#pragma once
#include "../SLISC/time.h"

#ifdef NDEBUG // release mode
#define SLS_TIME_H_ERR(str) SLS_ERR(str)
#else // debug mode
void SLS_TIME_H_ERR(const std::string &str) {}
#endif

// test time utilities
void test_time()
{
#ifndef NDEBUG
    std::cout << "test_time() : error not reported in debug mode!" << std::endl;
#endif
    using namespace slisc;
    Timer t; CPUTimer cput;
    // cpu time
    cput.tic();
    if (cput.toc() > 0.1)
        SLS_TIME_H_ERR("failed!");

    // natural time
    t.tic(); pause(0.114);
    if (abs(t.toc() - 0.114) > 1e-4) SLS_TIME_H_ERR("failed!");
}
