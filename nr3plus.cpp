# include "nr3plus.h"
using std::cout; using std::endl; using std::vector;

// display scalar, vector, or matrix
// not using template (or functions will not be available during debug)

// version 1

std::chrono::steady_clock::time_point tic_time_start; // for tic() toc()
vector<std::chrono::steady_clock::time_point> tic_time_starts(20);

void disp(VecUchar_I &v)
{
	Long i, n{ v.size() };
	Int precision = 4;
	auto oldPrecision = cout.precision();
	cout.precision(precision);
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(v) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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
	if (numel(a) == 0) cout << "empty";
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

void dft(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax)
{
	Long i, j, k, Nx = X.nrows(), Nc = X.ncols();
	Doub dk = (kmax-kmin)/(Nk-1), dx = (xmax-xmin)/(Nx-1);
	const Comp *pxi;
	Comp *pyj, factor, expo, dexpo;
	Y.resize(Nk, Nc); Y = 0.;
	for (j = 0; j < Nk; ++j) {
		pyj = Y[j];
		expo = exp(Comp(0, -(kmin + dk*j)*(xmin-dx)));
		dexpo = exp(Comp(0, -(kmin + dk*j)*dx));
		for (i = 0; i < Nx; ++i) {
			pxi = X[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pyj[k] += expo*pxi[k];
		}
	}
}

void idft(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax)
{
	Long i, j, k, Nk = Y.nrows(), Nc = Y.ncols();
	Doub dk = (kmax-kmin)/(Nk-1), dx = (xmax-xmin)/(Nx-1);
	const Comp *pyi;
	Comp *pxj, factor, expo, dexpo;
	X.resize(Nx, Nc); X = 0.;
	for (j = 0; j < Nx; ++j) {
		pxj = X[j];
		expo = exp(Comp(0, (xmin + dx*j)*(kmin-dk)));
		dexpo = exp(Comp(0, (xmin + dx*j)*dk));
		for (i = 0; i < Nk; ++i) {
			pyi = Y[i];
			expo *= dexpo;
			for (k = 0; k < Nc; ++k)
				pxj[k] += expo*pyi[k];
		}
	}
}