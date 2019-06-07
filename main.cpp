// comprehensive test of SLISC
#include "test/test_all.h"
// using namespace slisc;

#ifdef _MSC_VER
slisc::turn_on_floating_exceptions yes_turn_on_floating_exceptions;
#endif

int main()
{
	test_all();
}
