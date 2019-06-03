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

	// slice column from column major matrix
	{
		CmatInt a(3, 4);
		VecInt vc(3);
		SvecInt svc;
		for (Long j = 1; j < 4; ++j) {
			slice_col(svc, a, j);
			vc = svc;
			for (Long i = 1; i < 3; ++i) {
				if (svc[i] != a(i, j))
					SLS_ERR("failed!");
				if (vc[i] != svc[i])
					SLS_ERR("failed!");
			}
		}
	}
	// slice row from row major matrix
	{
		MatInt a(3, 4);
		VecInt vr(4);
		SvecInt svr;
		for (Long i = 1; i < 3; ++i) {
			slice_row(svr, a, i);
			vr = svr;
			for (Long j = 1; j < 3; ++j) {
				if (svr[j] != a(i, j))
					SLS_ERR("failed!");
				if (vr[j] != svr[j])
					SLS_ERR("failed!");
			}
		}
	}

	// slice row from column major matrix
	{
		CmatInt a(3, 4);
		linspace(a, 1, 12);
		VecInt v(4);
		for (Long i = 0; i < 3; ++i) {
			v = slice_row(a, i);
			for (Long j = 0; j < 4; ++j) {
				if (v[j] != a(i, j))
					SLS_ERR("failed!");
			}
			v *= 2;
			DvecInt slice; slice_row(slice, a, i);
			slice = v;
			for (Long j = 0; j < 4; ++j) {
				if (v[j] != a(i, j))
					SLS_ERR("failed!");
			}
		}
	}

	// slice a3(i,j,:)
	{
		Cmat3Int a3(2, 3, 4);
		linspace(a3, 1, 24);
		DvecInt slice;
		for (Long i = 0; i < 2; ++i) {
			for (Long j = 0; j < 3; ++j) {
				slice_dim3(slice, a3, i, j);
				for (Long k = 0; k < 4; ++k) {
					if (slice[k] != a3(i, j, k))
						SLS_ERR("failed!");
				}
			}
		}
	}
	{
		Cmat3Int a3(2, 2, 2);
		linspace(a3, 1, 8);
		DvecInt slice;
		slice_dim3(slice, a3, 1, 1);
		if (slice[0] != 4 || slice[1] != 8)
			SLS_ERR("failed!");
		slice /= 2;
		if (slice[0] != 2 || slice[1] != 4)
			SLS_ERR("failed!");
	}
}
