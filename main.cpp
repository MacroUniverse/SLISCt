// comprehensive test of SLISC
#include "test/test_all.h"
// using namespace slisc;

#ifdef _MSC_VER
slisc::turn_on_floating_exceptions yes_turn_on_floating_exceptions;
slisc::set_windows_console_utf8 yes_set_windows_console_utf8;
#endif

int main()
{
	test_all();
}
