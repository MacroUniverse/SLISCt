// display vectors and matrices
// header-only version, in namespace slisc
// unable to be called in debugger

#pragma once
#include "slisc.h"

namespace slisc {

// version 1
void disp(VecUchar_I &v);
void disp(VecInt_I &v);
void disp(VecDoub_I &v);
void disp(VecComp_I &v);
void disp(MatUchar_I &a);
void disp(MatInt_I &a);
void disp(MatDoub_I &a);
void disp(MatComp_I &a);
void disp(Mat3Doub_I &a);
void disp(Mat3Comp_I &a);
// version 2
void disp(VecUchar_I &v, Int_I precision);
void disp(VecInt_I &v, Int_I precision);
void disp(VecDoub_I &v, Int_I precision);
void disp(VecComp_I &v, Int_I precision);
void disp(MatUchar_I &a, Int_I precision);
void disp(MatInt_I &a, Int_I precision);
void disp(MatDoub_I &a, Int_I precision);
void disp(MatComp_I &a, Int_I precision);
void disp(Mat3Doub_I &a, Int_I precision);
void disp(Mat3Comp_I &a, Int_I precision);
// version 3
void disp(VecUchar_I &v, Long_I start, Long_I n);
void disp(VecInt_I &v, Long_I start, Long_I n);
void disp(VecDoub_I &v, Long_I start, Long_I n);
void disp(VecComp_I &v, Long_I start, Long_I n);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3);
// version 4
void disp(VecUchar_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecInt_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecDoub_I &v, Long_I start, Long_I n, Int_I precision);
void disp(VecComp_I &v, Long_I start, Long_I n, Int_I precision);
void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision);
void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);
void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision);

// implementation


// version 1

inline void disp(VecUchar_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << (Int)v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecInt_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecDoub_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecComp_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatUchar_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << (Int)a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatInt_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatDoub_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatComp_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Doub_I &a)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = 0; k < q; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Comp_I &a)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = 0; k < q; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

// version 2

inline void disp(VecUchar_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << (Int)v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecInt_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecDoub_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecComp_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatUchar_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << (Int)a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatInt_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatDoub_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatComp_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Doub_I &a, Int_I precision)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = 0; k < q; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Comp_I &a, Int_I precision)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = 0; k < q; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

// version 3

inline void disp(VecUchar_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << (Int)v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecInt_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecDoub_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecComp_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << (Int)a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{
	Long i, j, k;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = start3; k < start3+n3; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = start1; i < start1+n1; ++i) {
			for (j = start2; j < start2+n2; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{
	Long i, j, k;
	Int precision{ 4 };
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = start3; k < start3 + n3; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

// version 4

inline void disp(VecUchar_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << (Int)v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecInt_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecDoub_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(VecComp_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (v.size() == 0) std::cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		std::cout << v[i] << "   ";
	}
	std::cout << std::endl << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << (Int)a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			std::cout << a[i][j] << "   ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{
	Long i, j, k;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = start3; k < start3 + n3; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

inline void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{
	Long i, j, k;
	auto oldPrecision = std::cout.precision();
	std::cout.precision(precision);
	if (a.size() == 0) std::cout << "empty";
	else
	for (k = start3; k < start3 + n3; ++k) {
		std::cout << "(:, :, " << k << ")" << std::endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				std::cout << a[i][j][k] << "   ";
			}
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}
	std::cout.precision(oldPrecision);
}

} // namespace slisc
