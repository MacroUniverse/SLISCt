#include "../SLISC/sparse.h"

inline void test_sparse()
{
	using namespace slisc;
	// default constructor
	{
		McooDoub a;
		if (a.size() != 0 || a.nrows() != 0 || a.ncols() != 0 || a.nnz() != 0)
			error("failed!");
	}

	// 3 arg constructor
	{
		McooDoub a(3,4,5);
		if (a.nrows() != 3 || a.ncols() != 4 || a.size() != 5 || a.nnz() != 0)
			error("failed!");
	}

	// push, indexing, trim
	{
		McooDoub a(3, 4, 5);
		a.push(1., 0, 0);
		if (a.nnz() != 1) error("failed!");
		if (a(0) != 1.) error("failed!");
		if (a.row(0) != 0 || a.col(0) != 0) error("failed!");
		a.push(2., 1, 2);
		if (a.nnz() != 2) error("failed!");
		if (a(1) != 2.) error("failed!");
		if (a.row(1) != 1 || a.col(1) != 2) error("failed!");
		a(1) = 2.2; a.row(1) = 2; a.col(1) = 1;
		if (a.nnz() != 2) error("failed!");
		if (a(1) != 2.2) error("failed!");
		if (a.row(1) != 2 || a.col(1) != 1) error("failed!");
		a.push(3., 1, 3);

		a.trim(2);
		if (a.nnz() != 2) error("failed!");
		if (a(1) != 2.2) error("failed!");
		if (a.row(1) != 2 || a.col(1) != 1) error("failed!");
		a.trim(0);
		if (a.nnz() != 0) error("failed!");
	}
	
	{
		// copy assignment
		Long i, j, k = 0;
		McooDoub a(3, 4, 12), a1;
		for (i = 0; i < 3; ++i)
		for (j = 0; j < 4; ++j) {
			++k;
			a.push(k, i, j);
		}
		a1 = a;
		if (a1.nrows() != a.nrows() || a1.ncols() != a.ncols() ||
			a1.nnz() != a.nnz() || a1.size() != a.size())
			error("failed!");
		k = 0;
		for (i = 0; i < 3; ++i)
		for (j = 0; j < 4; ++j) {
			if (a1.row(k) != a.row(k)) error("failed!");
			if (a1.col(k) != a.col(k)) error("failed!");
			if (a1(k) != a(k)) error("failed!");
			++k;
		}

		// matrix - vector multiplication
		VecDoub v(3), v1(4);
		v1(0) = 1.; v1(1) = 2.; v1(2) = 3.; v1(3) = 4.;
		mul(v, a, v1);
		if (abs(v(0)-30.) + abs(v(1)-70.) + abs(v(2)-110.) > 1e-15)
			error("failed!");
	}
}
