# include "nr3plus.h"
using std::cout; using std::endl;

// display scalar, vector, or matrix
// not using template (or functions will not be available during debug)

// version 1

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

// version 2

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

// version 3
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

// version 4
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