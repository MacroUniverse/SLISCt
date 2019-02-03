#pragma once
#include "../SLISC/disp.h"
#include "../SLISC/sparse.h"

inline void test_sparse()
{
	using namespace slisc;
	// default constructor
	{
		McooDoub a;
		if (a.size() != 0 || a.nrows() != 0 || a.ncols() != 0 || a.nnz() != 0)
			error("failed!");
		McoohComp b;
		if (b.size() != 0 || b.nrows() != 0 || b.ncols() != 0 || b.nnz() != 0)
			error("failed!");
	}

	// 2 arg constructor
	{
		McooDoub a(3, 4);
		if (a.nrows() != 3 || a.ncols() != 4 || a.size() != 0 || a.nnz() != 0)
			error("failed!");
		McoohComp b(4, 4);
		if (b.nrows() != 4 || b.ncols() != 4 || b.size() != 0 || b.nnz() != 0)
			error("failed!");
	}

	// 3 arg constructor
	{
		McooDoub a(3, 4, 5);
		if (a.nrows() != 3 || a.ncols() != 4 || a.size() != 5 || a.nnz() != 0)
			error("failed!");
		McoohComp b(3, 3, 5);
		if (b.nrows() != 3 || b.ncols() != 3 || b.size() != 5 || b.nnz() != 0)
			error("failed!");
	}

	{
		// push
		McooDoub a(3, 4, 5); McoohComp b(4, 4, 5);
		a.push(1., 0, 0); b.push(Comp(1., 1.), 0, 0);
		if (a.nnz() != 1 || b.nnz() != 1) error("failed!");
		if (a(0) != 1. || b(0) != Comp(1.,1.)) error("failed!");
		if (a.row(0) != 0 || a.col(0) != 0) error("failed!");
		if (b.row(0) != 0 || b.col(0) != 0) error("failed!");
		a.push(2., 1, 2); b.push(Comp(2., 2.), 1, 2);
		if (a.nnz() != 2 || b.nnz() != 2) error("failed!");
		if (a(1) != 2. || b(1) != Comp(2.,2.)) error("failed!");
		if (a.row(1) != 1 || a.col(1) != 2) error("failed!");
		if (b.row(1) != 1 || b.col(1) != 2) error("failed!");
		a.push(3., 1, 3); b.push(Comp(3., 3.), 1, 3);

		// single indexing
		a(1) = 2.2; b(1) = Comp(2.2, 2.2);
		if (a(1) != 2.2 || b(1) != Comp(2.2,2.2)) error("failed!");

		// double indexing (element must exist)
		if (a(0, 0) != 1. || b(0, 0) != Comp(1.,1.)) error("failed!");
		if (a(1, 2) != 2.2 || b(1, 2) != Comp(2.2,2.2)) error("failed!");
		if (a(1, 3) != 3. || b(1, 3) != Comp(3.,3.)) error("failed!");
		a(1, 3) = 3.3;
		b.ref(1, 3) = Comp(3.3, 3.3);
		if (a(1, 3) != 3.3 || b(1, 3) != Comp(3.3, 3.3)) error("failed!");

		// set
		a.set(4.4, 1, 3); b.set(Comp(4.4, 4.4), 1, 3);
		if (a(1, 3) != 4.4 || b(1, 3) != Comp(4.4, 4.4)) error("failed!");
		a.set(3.14, 2, 0); b.set(Comp(3.14,-3.14), 0, 2);
		if (a(2, 0) != 3.14 || b(2, 0) != Comp(3.14, 3.14)) error("failed!");

		// trim
		a.trim(2); b.trim(2);
		if (a.nnz() != 2 || b.nnz() != 2) error("failed!");
		if (a(1) != 2.2 || b(1) != Comp(2.2,2.2)) error("failed!");
		if (a.row(1) != 1 || a.col(1) != 2) error("failed!");
		if (b.row(1) != 1 || b.col(1) != 2) error("failed!");
		a.trim(0); b.trim(0);
		if (a.nnz() != 0 || b.nnz() != 0) error("failed!");
	}
	
	{
		// copy assignment
		Long i, j, k = 0;
		McooDoub a(4, 4, 16), a1;
		McoohComp b(4, 4, 10), b1;
		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j) {
				++k;
				a.push(k, i, j);
				if (i >= j) b.push(Comp(k, k), i, j);
			}
		for (i = 0; i < b.size(); ++i) {
			if (b.row(i) > b.col(i))
				error("failed!");
		}
		a1 = a; b1 = b;
		if (a1.nrows() != a.nrows() || a1.ncols() != a.ncols() ||
			a1.nnz() != a.nnz() || a1.size() != a.size())
			error("failed!");
		if (b1.nrows() != b.nrows() || b1.ncols() != b.ncols() ||
			b1.nnz() != b.nnz() || b1.size() != b.size())
			error("failed!");
		for (i = 0; i < a.size(); ++i) {
			if (a1.row(i) != a.row(i)) error("failed!");
			if (a1.col(i) != a.col(i)) error("failed!");
			if (a1(i) != a(i)) error("failed!");
			++k;
		}
		for (i = 0; i < b.size(); ++i) {
			if (b1.row(i) != b.row(i)) error("failed!");
			if (b1.col(i) != b.col(i)) error("failed!");
			if (b1(i) != b(i)) error("failed!");
			++k;
		}

		// copy from sparse to dense matrix
		MatDoub c; MatComp d;
		c = a; d = b;
		CmatDoub e; CmatComp f;
		e = a; f = b;
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				if (c(i, j) != a(i, j) || e(i, j) != a(i, j)) {
					error("failed!");
				}
				if (d(i, j) != b(i, j) || f(i, j) != b(i, j)) {
					error("failed!");
				}
			}
		}

		// matrix - vector multiplication
		VecDoub y, y1, x(4);
		linspace(x, 1., 4.);
		mul(y, a, x); mul(y1, c, x);
		if (y != y1) error("failed!");

		VecComp v, v1;
		mul(v, b, x); mul(v1, d, x);
		if (v != v1) error("failed!");
	}
}