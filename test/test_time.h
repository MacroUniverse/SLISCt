#include "../SLISC/time.h"

#ifdef NDEBUG // release mode
#define time_h_error(str) error(str)
#else // debug mode
void time_h_error(const std::string &str) {}
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
	cput.tic(); pause(0.114);
	if (abs(cput.toc() - 0.114) > 2e-4) time_h_error("failed!");

	// natural time
	Int ind;
	t.tic(); pause(0.114);
	if (abs(t.toc() - 0.114) > 1e-4) time_h_error("failed!");
}
