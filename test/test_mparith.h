#include "../SLISC/mparith.h"

inline void test_mparith()
{
	using std::cout; using std::endl; //using std::conj;
	using namespace slisc;
	MParith a;
	std::string str = a.mppi(20);
	if (str != "3.1415926535897932384626433832795028841971693993")
		error("failed!");
}
