# include "nr3plus.h"
using std::cout; using std::endl;

// display scalar, vector, or matrix
// not using template (or functions will not be available during debug)

// version 1

void disp(VecUchar_I &v)
{
	int i, n{ v.size() }, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v)
{
	int i, n{ v.size() }, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v)
{
	int i, n{ v.size() }, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComplex_I &v)
{
	int i, n{ v.size() }, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComplex_I &a)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3DDoub_I &a)
{
	int i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = 0; k < q; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

void disp(Mat3DComplex_I &a)
{
	int i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() }, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = 0; k < q; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

// version 2

void disp(VecUchar_I &v, const int precision)
{
	int i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, const int precision)
{
	int i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, const int precision)
{
	int i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComplex_I &v, const int precision)
{
	int i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, const int precision)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, const int precision)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, const int precision)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComplex_I &a, const int precision)
{
	int i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3DDoub_I &a, const int precision)
{
	int i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = 0; k < q; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

void disp(Mat3DComplex_I &a, const int precision)
{
	int i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = 0; k < q; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = 0; i < m; ++i) {
			for (j = 0; j < n; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

// version 3

void disp(VecUchar_I &v, const int start, const int n)
{
	int i, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, const int start, const int n)
{
	int i, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, const int start, const int n)
{
	int i, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComplex_I &v, const int start, const int n)
{
	int i, precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, const int start1, const int start2, const int n1, const int n2)
{
	int i, j, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, const int start1, const int start2, const int n1, const int n2)
{
	int i, j, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, const int start1, const int start2, const int n1, const int n2)
{
	int i, j, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComplex_I &a, const int start1, const int start2, const int n1, const int n2)
{
	int i, j, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3DDoub_I &a, const int start1, const int start2, const int start3, const int n1, const int n2, const int n3)
{
	int i, j, k, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = start3; k < start3+n3; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = start1; i < start1+n1; ++i) {
			for (j = start2; j < start2+n2; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

void disp(Mat3DComplex_I &a, const int start1, const int start2, const int start3, const int n1, const int n2, const int n3)
{
	int i, j, k, precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = start3; k < start3 + n3; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

// version 4

void disp(VecUchar_I &v, const int start, const int n, const int precision)
{
	int i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, const int start, const int n, const int precision)
{
	int i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, const int start, const int n, const int precision)
{
	int i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComplex_I &v, const int start, const int n, const int precision)
{
	int i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, const int start1, const int start2, const int n1, const int n2, const int precision)
{
	int i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, const int start1, const int start2, const int n1, const int n2, const int precision)
{
	int i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, const int start1, const int start2, const int n1, const int n2, const int precision)
{
	int i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComplex_I &a, const int start1, const int start2, const int n1, const int n2, const int precision)
{
	int i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3DDoub_I &a, const int start1, const int start2, const int start3, const int n1, const int n2, const int n3, const int precision)
{
	int i, j, k;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = start3; k < start3 + n3; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}

void disp(Mat3DComplex_I &a, const int start1, const int start2, const int start3, const int n1, const int n2, const int n3, const int precision)
{
	int i, j, k;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	for (k = start3; k < start3 + n3; ++k) {
		cout << "(:, :, " << k << ")" << endl;
		for (i = start1; i < start1 + n1; ++i) {
			for (j = start2; j < start2 + n2; ++j) {
				cout << a[i][j][k] << "   ";
			}
			cout << endl;
		}
		cout << endl;
	}
	cout.precision(oldPrecision);
}
