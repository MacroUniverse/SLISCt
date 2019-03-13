#include "../SLISC/slice.h"
#include "../SLISC/disp.h"

void test_slice()
{
	using namespace slisc;
	// slice(v, ptr)
	CmatInt a(3,3);
	linspace(a, 1, 9);
	VecInt v(3);
	slice(v, a.ptr());
	if (v[0] != 1 || v[1] != 2 || v[2] != 3) error("failed!");
	v[0] = v[1] = v[2] = 0;
	if (v[0] != 0 || v[1] != 0 || v[2] != 0) error("failed!");
	slice(v, &a(0, 1), 2);
	if (v.size() != 2) error("failed!");
	if (v[0] != 4 || v[1] != 5) error("failed!");
	rm_slice(v);
}
