#include "../SLISC/slice_arith.h"
#include "../SLISC/disp.h"

void test_slice()
{
	using namespace slisc;
	// === test "Svector<>" ===
	CmatInt a(3,3);
	linspace(a, 1, 9);
	SvecInt v(a.ptr(), 3);
	if (v[0] != 1 || v[1] != 2 || v[2] != 3) SLS_ERR("failed!");
	v[0] = v[1] = v[2] = 0;
	if (v[0] != 0 || v[1] != 0 || v[2] != 0) SLS_ERR("failed!");
	v.set(&a(0, 1), 2);
	if (v.size() != 2) SLS_ERR("failed!");
	if (v[0] != 4 || v[1] != 5) SLS_ERR("failed!");
	v.set_ptr(&a(0, 2));
	if (v.size() != 2) SLS_ERR("failed!");
	if (v[0] != 7 || v[1] != 8 ) SLS_ERR("failed!");
	v.set_size(3);
	if (v.size() != 3) SLS_ERR("failed!");
	if (v[0] != 7 || v[1] != 8 || v[2] != 9) SLS_ERR("failed!");

	// arithmetics
	if (sum(v) != 24) SLS_ERR("failed!");
	v += 3;
	if (sum(v) !=33) SLS_ERR("failed!");
}
