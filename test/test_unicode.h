#pragma once
#include "../SLISC/unicode.h"

inline void test_unicode()
{
	using namespace slisc;
#ifdef _MSC_VER
	cout << "this will only display correctly in Visual C++ or on Linux OS :" << endl;
	cout << u8"显示一些中文" << endl;
#endif
}
