#include "../SLISC/time.h"

#ifdef NDEBUG // release mode
#define time_h_error(str) error(str)
#else // debug mode
#define time_h_error(str)
#endif

// test time utilities
void test_time()
{
#ifndef NDEBUG
	std::cout << "test_time() : error not reported in debug mode!" << std::endl;
#endif
	using namespace slisc;
	// cpu time
	ctic(); pause(0.114);
	if (abs(ctoc() - 0.114) > 2e-4) time_h_error("failed!");

	// physical time
	Int ind;
	tic(); pause(0.114);
	if (abs(toc() - 0.114) > 1e-4) time_h_error("failed!");
	tic(ind); pause(0.114);
	if (abs(toc(ind) - 0.114) > 1e-4) time_h_error("failed!");
	if (ind != 0) time_h_error("failed!");
	tic(ind); pause(0.114);
	if (abs(toc(ind) - 0.114) > 1e-4) time_h_error("failed!");
	if (ind != 1) time_h_error("failed!");
}
