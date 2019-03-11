// save vectors and matrices defined in "nr3.h" to ".mat" or ".matt" files.
// see README.txt for details
// class types: Doub=0, Comp=1, Int=2, Uchar=3.

#pragma once
//#define MATFILE_BINARY
//#define MATFILE_DUAL

#ifndef MATFILE_PRECISION
#define MATFILE_PRECISION 16
#endif

#include "../SLISC/slisc.h"

namespace slisc {

// matt file object
struct MATTFile;

// open a file (read: rw = 'r', write: rw = 'w')
MATTFile *mattOpen(Str fname, Char_I *rw);

// close a file;
void mattClose(MATTFile *pfile);

// save a scalar / matrix
// void mattsave(const T &s, Str_I varname, MATTFile *pfile);

// load a scalar / matrix
// void mattload(T &I, Str_I varname, MATTFile *pfile);

// ========== Implementation ============

// MATTFile class for text mode
struct MATTFile {
	char rw; // 'r' for read 'w' for write
	ifstream in; // read file
	ofstream out; // write file
	Int n; // variable numbers
	vector<Str> name; // variable names
	vector<Int> type; // variable types
	vector<vector<Long>> size; // variable dimensions
	vector<Long> ind; // variable positions (line indices)
};

// read the next variable after previous '\n'
Long scanInverse(std::ifstream &fin)
{
	Char c;
	Long ind, i, N;

	ind = fin.tellg();
	for (i = 2; i < 100; ++i) {
		fin.seekg(ind - i); c = fin.get();
		if (c == '\n') break;
	}
	fin >> N;
	fin.seekg(ind - i);
	return N;
}

// get var names and positions from the end of the file
// after return, pfile->ind[i] points to the first matrix element;
inline void getprofile(MATTFile *pfile)
{
	Int i, j, n, temp;
	std::vector<Long> size;
	Str name;
	std::ifstream &fin = pfile->in;

	// read number of variables and their positions
	fin.seekg(0, fin.end);
	pfile->n = (Int)scanInverse(fin);
	for (i = 0; i < pfile->n; ++i)
		pfile->ind.push_back(scanInverse(fin));

	// loop through each variable
	for (i = 0; i < pfile->n; ++i) {
		fin.seekg(pfile->ind[i]);
		// read var name
		fin >> n;
		name.resize(0);
		for (j = 0; j < n; ++j) {
			fin >> temp; name.push_back((char)temp);
		}
		pfile->name.push_back(name);
		// read var type
		fin >> temp; pfile->type.push_back(temp);
		// read var dim
		fin >> n;
		size.resize(0);
		for (j = 0; j < n; ++j) {
			fin >> temp; size.push_back(temp);
		}
		pfile->size.push_back(size);
		pfile->ind[i] = fin.tellg();
	}
}

MATTFile *mattOpen(Str fname, Char_I *rw)
{
	// must open file in binary mode, otherwise, '\n' will be written as "\r\n"
	// and seekg() will not work the same in linux.

	MATTFile* pfile = new MATTFile;
	if (rw[0] == 'w') {
		pfile->rw = 'w';
		pfile->n = 0;
		pfile->out = std::ofstream(fname, std::ios_base::binary);
		pfile->out.precision(MATFILE_PRECISION);
	}
	else {
		pfile->rw = 'r';
		pfile->in = std::ifstream(fname, std::ios_base::binary);
		if (!pfile->in)
			error("error: file not found: ");
		pfile->in.precision(17);
		getprofile(pfile); // get var names
	}
	return pfile;
}

inline void mattClose(MATTFile* pfile)
{
	Llong i;
	if (pfile->rw == 'w') {
		std::ofstream &fout = pfile->out;
		// write position of variables
		for (i = pfile->ind.size() - 1; i >= 0; --i)
			fout << pfile->ind[i] << "\n";
		// write number of variables
		fout << pfile->n;
		pfile->out.close();
	}
	else {
		pfile->in.close();
	}
	delete pfile;
}

template <class T, SLS_IF(
	is_Uchar<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
void mattsave(const T &s, Str_I varname, MATTFile *pfile)
{
	Long i, n;
	ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << '\n';
	}
	// write data type info
	if (is_Uchar<T>()) fout << 3 << '\n';
	else if (is_Int<T>()) fout << 2 << '\n';
	else if (is_Doub<T>()) fout << 0 << '\n';
	else if (is_Comp<T>()) fout << 1 << '\n';
	else error("unhandled case!");
	// write dimension info
	fout << 0 << '\n';
	// write matrix data
	if (!is_comp<T>()) {
		fout << to_num(s) << '\n';
	}
	else if (is_Comp<T>()) {
		if (imag(s) == 0)
			fout << real(s) << '\n';
		else if (imag(s) < 0)
			fout << real(s) << imag(s) << "i\n";
		else
			fout << real(s) << '+' << imag(s) << "i\n";
	}
	else
		error("unhandled case!");
}

template <class T, SLS_IF(
	is_Int<T>() || is_Uchar<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void mattsave(const Vector<T> &v, Str_I varname, MATTFile *pfile)
{
	Long i, n;
	ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << '\n';
	}
	// write data type info
	if (is_Uchar<T>()) fout << 3 << '\n';
	else if (is_Int<T>()) fout << 2 << '\n';
	else if (is_Doub<T>()) fout << 0 << '\n';
	else if (is_Comp<T>()) fout << 1 << '\n';
	else error("unhandled");
	// write dimension info
	n = v.size();
	fout << 1 << '\n' << n << '\n';
	// write matrix data
	for (i = 0; i < n; ++i) {
		if (!is_comp<T>())
			fout << to_num(v[i]) << '\n';
		else {
			auto cr = real(v[i]), ci = imag(v[i]);
			if (ci == 0)
				fout << cr << '\n';
			else if (ci < 0)
				fout << cr << ci << "i\n";
			else
				fout << cr << '+' << ci << "i\n";
		}
	}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	(is_Matrix<Tm>() || is_Cmat<Tm>()) &&
	(is_Uchar<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void mattsave(const Tm &a, Str_I varname, MATTFile *pfile,
	Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, m, n;
	ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << '\n';
	}
	// write data type info
	if (is_Uchar<T>()) fout << 3 << '\n';
	else if (is_Int<T>()) fout << 2 << '\n';
	else if (is_Doub<T>()) fout << 0 << '\n';
	else if (is_Comp<T>()) fout << 1 << '\n';
	else error("unhandled!");
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << '\n' << m << '\n' << n << '\n';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			if (!is_comp<T>())
				fout << to_num(a(step1*i, step2*j)) << '\n';
			else {
				auto c = a(step1*i, step2*j); auto cr = real(c), ci = imag(c);
				if (ci == 0)
					fout << cr << '\n';
				else if (ci < 0)
					fout << cr << ci << "i\n";
				else
					fout << cr << '+' << ci << "i\n";
			}
		}
}

template <class T, SLS_IF(
	is_Doub<T>() || is_Comp<T>()
)>
inline void mattsave(const Mat3d<T> &a, Str_I varname, MATTFile *pfile,
	Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1)
{
	Long i, j, k, m, n, q;
	ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << '\n';
	}
	// write data type info
	if (is_Doub<T>()) fout << 0 << '\n';
	else if (is_Comp<T>()) fout << 1 << '\n';
	// write dimension info
	m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
	q = (a.dim3() + step3 - 1) / step3;
	fout << 3 << '\n' << m << '\n' << n << '\n' << q << '\n';
	// write matrix data
	for (k = 0; k < q; ++k)
	for (j = 0; j < n; ++j)
	for (i = 0; i < m; ++i) {
		if constexpr (!is_comp<T>()) {
			fout << a(step1*i, step2*j, step3*k) << '\n';
		}
		else {
			auto c = a(step1*i, step2*j, step3*k); auto cr = real(c), ci = imag(c);
			if (ci == 0)
				fout << cr << '\n';
			else if (ci < 0)
				fout << cr << ci << "i\n";
			else
				fout << cr << '+' << ci << "i\n";
		}
	}
}

inline void mattsave(Mat3Doub_I &a, Str_I varname, MATTFile *pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 0 << '\n';
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
					fout << a(ind, step1*j, step2*k) << '\n';
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k)
					fout << a(step2*i, ind, step1*k) << '\n';
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i)
					fout << a(step1*i, step2*j, ind) << '\n';
		}
	}
	else
		error("illegal value of xyz");
}

inline void mattsave(Mat3Comp_I &a, Str_I varname, MATTFile *pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	Comp c; Doub cr, ci;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 1 << '\n';
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
				{
					c = a(ind, step1*j, step2*k); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << '\n';
					else if (ci < 0)
						fout << cr << ci << "i\n";
					else
						fout << cr << '+' << ci << "i\n";
				}
					
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k) {
					c = a(step2*i, ind, step1*k); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << '\n';
					else if (ci < 0)
						fout << cr << ci << "i\n";
					else
						fout << cr << '+' << ci << "i\n";
				}
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << '\n' << m << '\n' << n << '\n' << Nslice << '\n';
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i) {
					c = a(step1*i, step2*j, ind); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << '\n';
					else if (ci < 0)
						fout << cr << ci << "i\n";
					else
						fout << cr << '+' << ci << "i\n";
				}
		}
	}
	else
		error("illegal value of xyz");
}

// search variable in file by name
inline Int nameSearch(Str_I name, MATTFile *pfile)
{
	for (Int i = 0; i < pfile->n; ++i)
		if (name == pfile->name[i])
			return i;
	error("variable name not found!");
	return -1;
}

inline void scanComplex(Comp &c, std::ifstream &fin)
{
	Doub cr = 0, ci = 0;
	Uchar ch;
	fin >> cr;
	ch = fin.get();
	if (ch == '\n') {
		c = cr; return;
	}
	fin >> ci;
	if (ch == '-')
		ci *= -1.;
	c = Comp(cr, ci);
	fin.ignore(100, '\n');
}

template <class T, SLS_IF(
	is_Uchar<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void mattload(T &s, Str_I varname, MATTFile *pfile)
{
	Int i;
	ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim and var data
	if constexpr (is_Uchar<T>()) {
		if (pfile->type[i] != 3 || pfile->size[i].size() != 0)
			error("wrong type or dim!");
		Int temp; fin >> temp; s = (T)temp;
	}
	else if constexpr (is_Int<T>()) {
		// read var type and dim
		if (pfile->type[i] < 2 || pfile->size[i].size() != 0)
			error("wrong type or dim!");
		fin >> s;
	}
	else if constexpr (is_Doub<T>()) {
		if (pfile->type[i] == 1 || pfile->size[i].size() != 0)
			error("wrong type or dim!");
		fin >> s;
	}
	else if constexpr (is_Comp<T>()) {
		if (pfile->size[i].size() != 0)
			error("wrong type or dim!");
		scanComplex(s, fin);
	}
	else error("unhandled!");
}

template <class T, SLS_IF(
	is_Uchar<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void mattload(Vector<T> &v, Str_I varname, MATTFile *pfile)
{
	Long i, n, dim;
	ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (is_Uchar<T>()) {
		if (pfile->type[i] != 3 || dim != 1) error("wrong type or dim!");
	}
	else if (is_Int<T>()) {
		if (pfile->type[i] < 2 || dim != 1) error("wrong type or dim!");
	}
		
	else if (is_Doub<T>()) {
		if (pfile->type[i] == 1 || dim != 1) error("wrong type or dim!");
	}
		
	else if (is_Comp<T>()) {
		if (dim != 1) error("wrong type or dim!");
	}
	else error("unhandled!");

	n = pfile->size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i) {
		if constexpr (is_Uchar<T>()) {
			Int temp; fin >> temp;  v[i] = (Uchar)temp;
		}
		else if constexpr (is_Int<T>() || is_Doub<T>())
			fin >> v[i];
	}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	is_Matrix<Tm>() && (is_Uchar<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void mattload(Tm &a, Str_I varname, MATTFile *pfile)
{
	Long i, j, dim, m, n;
	ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (is_Uchar<T>()) {
		if (pfile->type[i] != 3 || dim != 2) error("wrong type or dim!");
	}
	else if (is_Int<T>()) {
		if (pfile->type[i] < 2 || dim != 2) error("wrong type or dim!");
	}
	else if (is_Doub<T>()) {
		if (pfile->type[i] == 1 || dim != 2) error("wrong type or dim!");
	}
	else if (is_Comp<T>()) {
		if (dim != 2) error("wrong type or dim!");
	}
	else error("unhandled!");

	m = pfile->size[i][0]; n = pfile->size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			if constexpr (is_Uchar<T>()) {
				Int temp; fin >> temp;  a(i, j) = (Uchar)temp;
			}
			else if constexpr (is_Int<T>() || is_Doub<T>())
				fin >> a(i, j);
			else if constexpr (is_comp<T>())
				scanComplex(a(i, j), fin);
			else error("unhandled!");
		}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	is_Mat3d<Tm>() && (is_Doub<T>() || is_Comp<T>())
)>
inline void mattload(Tm &a, Str_I varname, MATTFile *pfile)
{
	Long i, j, k, dim, m, n, q;
	ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (is_Doub<T>()) {
		if (pfile->type[i] == 1 || dim != 3) error("wrong type or dim!");
	}
	else if (is_Comp<T>()) {
		if (dim != 3) error("wrong type or dim!");
	}
	
	m = pfile->size[i][0]; n = pfile->size[i][1]; q = pfile->size[i][2];
	a.resize(m, n, q);
	// read var data
	for (k = 0; k < q; ++k)
		for (j = 0; j < n; ++j)
			for (i = 0; i < m; ++i) {
				if constexpr (is_Doub<T>())
					fin >> a(i, j, k);
				else if constexpr (is_Comp<T>())
					scanComplex(a(i, j, k), fin);
				else error("unhandled!");
			}
				
}

} // namespace slisc
