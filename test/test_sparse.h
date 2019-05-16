#pragma once
#include "../SLISC/arithmetic.h"
#include "../SLISC/sparse_arith.h"
#include "../SLISC/disp.h"

inline void test_sparse()
{
	using namespace slisc;
	// default constructor
	{
		McooDoub a(0,0);
		if (a.nnz() != 0 || a.nrows() != 0 || a.ncols() != 0 || a.nnz() != 0)
			SLS_ERR("failed!");
		McoohComp b(0,0);
		if (b.nnz() != 0 || b.nrows() != 0 || b.ncols() != 0 || b.nnz() != 0)
			SLS_ERR("failed!");
		DiagInt c(0);
		if (c.nnz() != 0 || c.nrows() != 0 || c.ncols() != 0)
			SLS_ERR("failed!");
	}

	// 1 arg constructor
	{
		DiagInt c(4);
		if (c.nnz() != 4 || c.nrows() != 4 || c.ncols() != 4)
			SLS_ERR("failed!");
	}

	// 2 arg constructor
	{
		McooDoub a(3, 4);
		if (a.nrows() != 3 || a.ncols() != 4 || a.capacity() != 0 || a.nnz() != 0)
			SLS_ERR("failed!");
		McoohComp b(4, 4);
		if (b.nrows() != 4 || b.ncols() != 4 || b.capacity() != 0 || b.nnz() != 0)
			SLS_ERR("failed!");
		DiagComp c(4, 4);
		if (c.nrows() != 4 || c.ncols() != 4 || c.nnz() != 4)
			SLS_ERR("failed!");
		DiagComp d(4, Comp(1,-3));
		if (d.nrows() != 4 || d.ncols() != 4 || d.nnz() != 4)
			SLS_ERR("failed!");
		if ((VecComp&)d != Comp(1, -3))
			SLS_ERR("failed!");
	}

	// 3 arg constructor
	{
		McooDoub a(3, 4, 5);
		if (a.nrows() != 3 || a.ncols() != 4 || a.capacity() != 5 || a.nnz() != 0)
			SLS_ERR("failed!");
		McoohComp b(3, 3, 5);
		if (b.nrows() != 3 || b.ncols() != 3 || b.capacity() != 5 || b.nnz() != 0)
			SLS_ERR("failed!");
		DiagInt c(5, 5, 1);
		if (c.nrows() != 5 || c.ncols() != 5 || c.nnz() != 5)
			SLS_ERR("failed!");
		if ((VecInt&)c != 1)
			SLS_ERR("failed!");
	}

	{
		// push
		McooDoub a(3, 4, 5); McoohComp b(4, 4, 5);
		a.push(1., 0, 0); b.push(Comp(1., 1.), 0, 0);
		if (a.nnz() != 1 || b.nnz() != 1) SLS_ERR("failed!");
		if (a(0) != 1. || b(0) != Comp(1.,1.)) SLS_ERR("failed!");
		if (a.row(0) != 0 || a.col(0) != 0) SLS_ERR("failed!");
		if (b.row(0) != 0 || b.col(0) != 0) SLS_ERR("failed!");
		a.push(2., 1, 2); b.push(Comp(2., 2.), 1, 2);
		if (a.nnz() != 2 || b.nnz() != 2) SLS_ERR("failed!");
		if (a(1) != 2. || b(1) != Comp(2.,2.)) SLS_ERR("failed!");
		if (a.row(1) != 1 || a.col(1) != 2) SLS_ERR("failed!");
		if (b.row(1) != 1 || b.col(1) != 2) SLS_ERR("failed!");
		a.push(3., 1, 3); b.push(Comp(3., 3.), 1, 3);

		// single indexing
		a(1) = 2.2; b(1) = Comp(2.2, 2.2);
		if (a(1) != 2.2 || b(1) != Comp(2.2,2.2)) SLS_ERR("failed!");

		// double indexing (read-only)
		if (a(0, 0) != 1. || b(0, 0) != Comp(1.,1.)) SLS_ERR("failed!");
		if (a(1, 2) != 2.2 || b(1, 2) != Comp(2.2,2.2)) SLS_ERR("failed!");
		if (a(1, 3) != 3. || b(1, 3) != Comp(3.,3.)) SLS_ERR("failed!");
		if (a(2, 3) != 0. || b(2, 3) != 0.) SLS_ERR("failed!");

		// ref()
		a.ref(1, 3) = 3.3;
		b.ref(1, 3) = Comp(3.3, 3.3);
		if (a(1, 3) != 3.3 || b(1, 3) != Comp(3.3, 3.3)) SLS_ERR("failed!");

		// set
		a.set(4.4, 1, 3); b.set(Comp(4.4, 4.4), 1, 3);
		if (a(1, 3) != 4.4 || b(1, 3) != Comp(4.4, 4.4)) SLS_ERR("failed!");
		a.set(3.14, 2, 0); b.set(Comp(3.14,-3.14), 0, 2);
		if (a(2, 0) != 3.14 || b(2, 0) != Comp(3.14, 3.14)) SLS_ERR("failed!");

		// trim
		a.trim(2); b.trim(2);
		if (a.nnz() != 2 || b.nnz() != 2) SLS_ERR("failed!");
		if (a(1) != 2.2 || b(1) != Comp(2.2,2.2)) SLS_ERR("failed!");
		if (a.row(1) != 1 || a.col(1) != 2) SLS_ERR("failed!");
		if (b.row(1) != 1 || b.col(1) != 2) SLS_ERR("failed!");
		a.trim(0); b.trim(0);
		if (a.nnz() != 0 || b.nnz() != 0) SLS_ERR("failed!");
	}
	
	// TODO: Diag

	{
		// copy assignment
		Long i, j, k = 0;
		McooDoub a(4, 4, 16), a1(0,0);
		McoohComp b(4, 4, 10), b1(0,0);
		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j) {
				++k;
				a.push(k, i, j);
				if (i >= j) b.push(Comp(k, k), i, j);
			}
		for (i = 0; i < b.nnz(); ++i) {
			if (b.row(i) > b.col(i))
				SLS_ERR("failed!");
		}
		a1 = a; b1 = b;
		if (a1.nrows() != a.nrows() || a1.ncols() != a.ncols() ||
			a1.nnz() != a.nnz() || a1.capacity() != a.capacity())
			SLS_ERR("failed!");
		if (b1.nrows() != b.nrows() || b1.ncols() != b.ncols() ||
			b1.nnz() != b.nnz() || b1.capacity() != b.capacity())
			SLS_ERR("failed!");
		for (i = 0; i < a.nnz(); ++i) {
			if (a1.row(i) != a.row(i)) SLS_ERR("failed!");
			if (a1.col(i) != a.col(i)) SLS_ERR("failed!");
			if (a1(i) != a(i)) SLS_ERR("failed!");
			++k;
		}
		for (i = 0; i < b.nnz(); ++i) {
			if (b1.row(i) != b.row(i)) SLS_ERR("failed!");
			if (b1.col(i) != b.col(i)) SLS_ERR("failed!");
			if (b1(i) != b(i)) SLS_ERR("failed!");
			++k;
		}

		// copy from sparse to dense matrix
		MatDoub c(a.nrows(), a.ncols()); MatComp d(b.nrows(), b.ncols());
		c = a; d = b;
		CmatDoub e(a.nrows(), a.ncols()); CmatComp f(b.nrows(),b.ncols());
		e = a; f = b;
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				if (c(i, j) != a(i, j) || e(i, j) != a(i, j)) {
					SLS_ERR("failed!");
				}
				if (d(i, j) != b(i, j) || f(i, j) != b(i, j)) {
					SLS_ERR("failed!");
				}
			}
		}
	}

	// matrix - vector multiplication
	{
		// for McooDoub
		{
			McooDoub a(4, 4, 16);
			for (Int i = 0; i < 16; ++i) {
				a.push(i, i % 4, i / 4);
			}
			CmatDoub a1(a.nrows(), a.ncols()); a1 = a;
			VecDoub v(a.nrows()), v1(a1.nrows()), x(4);
			linspace(x, 1., 4.);
			mul(v, a, x); mul(v1, a1, x);
			if (v != v1) SLS_ERR("failed!");
		}
		
		
		// for McooComp
		{
			McooComp b(4, 4, 16);
			for (Int i = 0; i < 16; ++i) {
				b.push(Comp(i, i + 1), i % 4, i / 4);
			}
			MatComp b1(b.nrows(), b.ncols()); b1 = b;
			VecComp v(b.nrows()), v1(b1.nrows()), x(4);
			linspace(x, Comp(1.,-1.), Comp(4.,-4.));
			mul(v, b, x); mul(v1, b1, x);
			if (v != v1) SLS_ERR("failed!");
		}
		
		// for McooImag
		{
			McooImag c(4, 4, 16);
			for (Int i = 0; i < 16; ++i) {
				c.push(Imag(i), i % 4, i / 4);
			}
			CmatImag c1(c.nrows(), c.ncols()); c1 = c;
			VecComp v(c.nrows()), v1(c1.nrows()), x(4);
			linspace(x, Comp(1.,-1.), Comp(4.,-4.));
			mul(v, c, x); mul(v1, c1, x);
			if (v != v1) SLS_ERR("failed!");
		}
	}

	// inf_norm
	{
		McooDoub a(3, 3, 10);
		a.push(1., 0, 0); a.push(-2., 1, 1); a.push(3., 2, 2);
		a.push(4., 0, 1); a.push(4., 1, 0);
		a.push(5., 1, 2); a.push(-5., 2, 1);
		if (norm_inf(a) != 11.) SLS_ERR("failed!");

		McoohComp b(3, 3, 10);
		b.push(1., 0, 0); b.push(-2., 1, 1); b.push(3., 2, 2);
		b.push(Comp(4.,5.), 0, 1); b.push(Comp(3.,2.), 1, 2);
		if (abs(norm_inf(b) - 12.0086755) > 1e-6) SLS_ERR("failed!");
	}

	// mul(cmat, cmat, diag)
	{
		CmatInt a(4, 5, 1);
		VecInt b(5); linspace(b, 1, 5);
		CmatInt c(4, 5);
		mul(c, a, (DiagInt)b);
		if (!shape_cmp(c, a)) SLS_ERR("failed!");
		if (!equals_to_vs(&c(0,0), 1, c.nrows())) SLS_ERR("failed!");
		if (!equals_to_vs(&c(0,1), 2, c.nrows())) SLS_ERR("failed!");
		if (!equals_to_vs(&c(0,2), 3, c.nrows())) SLS_ERR("failed!");
		if (!equals_to_vs(&c(0,3), 4, c.nrows())) SLS_ERR("failed!");
		if (!equals_to_vs(&c(0,4), 5, c.nrows())) SLS_ERR("failed!");
	}

	// mul(cmat, diag, cmat)
	{
		CmatDoub a(3, 3), b(3,3); linspace(a, 1, 9);
		VecDoub v(3); linspace(v, 1, 3);
		mul(b, diag(v), a);
		if (b[0] != 1 || b[1] != 4 || b[2] != 9
			|| b[3] != 4 || b[4] != 10 || b[5] != 18 ||
			b[6] != 7 || b[7] != 16 || b[8] != 27)
			SLS_ERR("failed!");
	}
}
