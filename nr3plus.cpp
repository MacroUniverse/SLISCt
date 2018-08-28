# include "nr3plus.h"
using std::cout; using std::endl; using std::vector;

// time utilities
std::chrono::steady_clock::time_point tic_time_start; // for tic() toc()
vector<std::chrono::steady_clock::time_point> tic_time_starts(20);
Llong ctic_time_start; // for ctic() ctoc()
std::vector<Llong> ctic_time_starts(20);

// display scalar, vector, or matrix
// not using template (or functions will not be available during debug)

// version 1

void disp(VecUchar_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComp_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComp_I &a)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3Doub_I &a)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(Mat3Comp_I &a)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(VecUchar_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComp_I &v, Int_I precision)
{
	Long i, n{ v.size() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = 0; i < n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComp_I &a, Int_I precision)
{
	Long i, j, m{ a.nrows() }, n{ a.ncols() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3Doub_I &a, Int_I precision)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(Mat3Comp_I &a, Int_I precision)
{
	Long i, j, k, m{ a.dim1() }, n{ a.dim2() }, q{ a.dim3() };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(VecUchar_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComp_I &v, Long_I start, Long_I n)
{
	Long i;
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2)
{
	Long i, j;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{
	Long i, j, k;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3)
{
	Long i, j, k;
	Int precision{ 4 };
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(VecUchar_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << (Int)v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecInt_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecDoub_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(VecComp_I &v, Long_I start, Long_I n, Int_I precision)
{
	Long i;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (v.size() == 0) cout << "empty";
	else
	for (i = start; i < start + n; ++i) {
		cout << v[i] << "   ";
	}
	cout << endl << endl;
	cout.precision(oldPrecision);
}

void disp(MatUchar_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << (Int)a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatInt_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatDoub_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(MatComp_I &a, Long_I start1, Long_I start2, Long_I n1, Long_I n2, Int_I precision)
{
	Long i, j;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
	for (i = start1; i < start1 + n1; ++i) {
		for (j = start2; j < start2 + n2; ++j) {
			cout << a[i][j] << "   ";
		}
		cout << endl;
	}
	cout << endl;
	cout.precision(oldPrecision);
}

void disp(Mat3Doub_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{
	Long i, j, k;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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

void disp(Mat3Comp_I &a, Long_I start1, Long_I start2, Long_I start3, Long_I n1, Long_I n2, Long_I n3, Int_I precision)
{
	Long i, j, k;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (a.size() == 0) cout << "empty";
	else
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
