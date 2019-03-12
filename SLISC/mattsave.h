// save vectors and matrices defined in "nr3.h" to ".mat" or ".matt" files.
// see README.txt for details
// class types: Doub=0, Comp=1, Int=2, Char=3.

#pragma once
//#define MATFILE_BINARY
//#define MATFILE_DUAL

#include "../SLISC/slisc.h"

namespace slisc {

// matt file object
struct Matt;

// ========== Implementation ============

// Matt class for text mode
class Matt {
public:
	char m_rw; // 'r' for read 'w' for write
	ifstream m_in; // read file
	ofstream m_out; // write file
	Int m_n; // variable numbers
	vector<Str> m_name; // variable names
	vector<Int> m_type; // variable types
	vector<vector<Long>> m_size; // variable dimensions
	vector<Long> m_ind; // variable positions (line indices)

	void get_profile();
	void open(Str fname, Char_I *rw, Int_I precision = 17);
	void close();
};

typedef const Matt &Matt_I;
typedef Matt &Matt_O, &Matt_IO;

// read the next variable after previous ' '
Long scanInverse(std::ifstream &fin)
{
	Char c;
	Long ind, i, N;

	ind = fin.tellg();
	for (i = 2; i < 100; ++i) {
		fin.seekg(ind - i); c = fin.get();
		if (c == ' ') break;
	}
	fin >> N;
	fin.seekg(ind - i);
	return N;
}

// get var names and positions from the end of the file
// after return, pfile.m_ind[i] points to the first matrix element;
inline void Matt::get_profile()
{
	Int i, j, n, temp;
	vector<Long> size;
	Str name;
	ifstream &fin = m_in;

	// read number of variables and their positions
	fin.seekg(0, fin.end);
	m_n = (Int)scanInverse(fin);
	for (i = 0; i < m_n; ++i)
		m_ind.push_back(scanInverse(fin));

	// loop through each variable
	for (i = 0; i < m_n; ++i) {
		fin.seekg(m_ind[i]);
		// read var name
		fin >> n;
		name.resize(0);
		for (j = 0; j < n; ++j) {
			fin >> temp; name.push_back((Char)temp);
		}
		m_name.push_back(name);
		// read var type
		fin >> temp; m_type.push_back(temp);
		// read var dim
		fin >> n;
		size.resize(0);
		for (j = 0; j < n; ++j) {
			fin >> temp; size.push_back(temp);
		}
		m_size.push_back(size);
		m_ind[i] = fin.tellg();
	}
}

void Matt::open(Str fname, Char_I *rw, Int_I precision)
{
	if (rw[0] == 'w') {
		m_rw = 'w';
		m_n = 0;
		m_out = std::ofstream(fname);
		m_out.precision(precision);
	}
	else {
		m_rw = 'r';
		m_in = std::ifstream(fname);
		if (!m_in)
			error("error: file not found: ");
		m_in.precision(17);
		get_profile(); // get var names
	}
}

inline void Matt::close()
{
	Llong i;
	if (m_rw == 'w') {
		std::ofstream &fout = m_out;
		// write position of variables
		for (i = m_ind.size() - 1; i >= 0; --i)
			fout << m_ind[i] << ' ';
		// write number of variables
		fout << m_n;
		m_out.close();
	}
	else {
		m_in.close();
	}
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
void save(const T &s, Str_I varname, Matt_IO pfile)
{
	Long i, n;
	ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << ' ';
	}
	// write data type info
	fout << type_num<T>() << ' ';
	// write dimension info
	fout << 0 << ' ';
	// write matrix data
	if (!is_comp<T>()) {
		fout << to_num(s) << ' ';
	}
	else if (is_Comp<T>()) {
		if (imag(s) == 0)
			fout << real(s) << ' ';
		else if (imag(s) < 0)
			fout << real(s) << imag(s) << "i ";
		else
			fout << real(s) << '+' << imag(s) << "i ";
	}
	else
		error("unhandled case!");
}

template <class T, SLS_IF(
	is_Int<T>() || is_Char<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void save(const Vector<T> &v, Str_I varname, Matt_IO pfile)
{
	Long i, n;
	ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << ' ';
	}
	// write data type info
	fout << type_num<T>() << ' ';
	// write dimension info
	n = v.size();
	fout << 1 << ' ' << n << ' ';
	// write matrix data
	for (i = 0; i < n; ++i) {
		if (!is_comp<T>())
			fout << to_num(v[i]) << ' ';
		else {
			auto cr = real(v[i]), ci = imag(v[i]);
			if (ci == 0)
				fout << cr << ' ';
			else if (ci < 0)
				fout << cr << ci << "i ";
			else
				fout << cr << '+' << ci << "i ";
		}
	}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	(is_Matrix<Tm>() || is_Cmat<Tm>()) &&
	(is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void save(const Tm &a, Str_I varname, Matt_IO pfile,
	Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, m, n;
	ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << ' ';
	}
	// write data type info
	fout << type_num<T>() << ' ';
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << ' ' << m << ' ' << n << ' ';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			if (!is_comp<T>())
				fout << to_num(a(step1*i, step2*j)) << ' ';
			else {
				auto c = a(step1*i, step2*j); auto cr = real(c), ci = imag(c);
				if (ci == 0)
					fout << cr << ' ';
				else if (ci < 0)
					fout << cr << ci << "i ";
				else
					fout << cr << '+' << ci << "i ";
			}
		}
}

template <class T, SLS_IF(
	is_Doub<T>() || is_Comp<T>()
)>
inline void save(const Mat3d<T> &a, Str_I varname, Matt_IO pfile,
	Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1)
{
	Long i, j, k, m, n, q;
	ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << ' ';
	}
	// write data type info
	fout << type_num<T>() << ' ';
	// write dimension info
	m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
	q = (a.dim3() + step3 - 1) / step3;
	fout << 3 << ' ' << m << ' ' << n << ' ' << q << ' ';
	// write matrix data
	for (k = 0; k < q; ++k)
	for (j = 0; j < n; ++j)
	for (i = 0; i < m; ++i) {
		if constexpr (!is_comp<T>()) {
			fout << a(step1*i, step2*j, step3*k) << ' ';
		}
		else {
			auto c = a(step1*i, step2*j, step3*k); auto cr = real(c), ci = imag(c);
			if (ci == 0)
				fout << cr << ' ';
			else if (ci < 0)
				fout << cr << ci << "i ";
			else
				fout << cr << '+' << ci << "i ";
		}
	}
}

inline void save(Mat3Doub_I &a, Str_I varname, Matt_IO pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	std::ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << ' ';
	}
	// write data type info
	fout << type_num<Doub>() << ' ';
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
					fout << a(ind, step1*j, step2*k) << ' ';
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k)
					fout << a(step2*i, ind, step1*k) << ' ';
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i)
					fout << a(step1*i, step2*j, ind) << ' ';
		}
	}
	else
		error("illegal value of xyz");
}

inline void save(Mat3Comp_I &a, Str_I varname, Matt_IO pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	Comp c; Doub cr, ci;
	std::ofstream &fout = pfile.m_out;
	++pfile.m_n; pfile.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << ' ';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << ' ';
	}
	// write data type info
	fout << type_num<Comp>() << ' ';
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
				{
					c = a(ind, step1*j, step2*k); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << ' ';
					else if (ci < 0)
						fout << cr << ci << "i ";
					else
						fout << cr << '+' << ci << "i ";
				}
					
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k) {
					c = a(step2*i, ind, step1*k); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << ' ';
					else if (ci < 0)
						fout << cr << ci << "i ";
					else
						fout << cr << '+' << ci << "i ";
				}
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << ' ' << m << ' ' << n << ' ' << Nslice << ' ';
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i) {
					c = a(step1*i, step2*j, ind); cr = real(c); ci = imag(c);
					if (ci == 0)
						fout << cr << ' ';
					else if (ci < 0)
						fout << cr << ci << "i ";
					else
						fout << cr << '+' << ci << "i ";
				}
		}
	}
	else
		error("illegal value of xyz");
}

// search variable in file by name
inline Int nameSearch(Str_I name, Matt_IO pfile)
{
	for (Int i = 0; i < pfile.m_n; ++i)
		if (name == pfile.m_name[i])
			return i;
	error("variable name not found!");
	return -1;
}

inline void scanComplex(Comp &c, std::ifstream &fin)
{
	Doub cr = 0, ci = 0;
	Char ch;
	fin >> cr;
	ch = fin.get();
	if (ch == ' ') {
		c = cr; return;
	}
	fin >> ci;
	if (ch == '-')
		ci *= -1.;
	c = Comp(cr, ci);
	fin.ignore(100, ' ');
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void load(T &s, Str_I varname, Matt_IO pfile)
{
	Int i;
	ifstream &fin = pfile.m_in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile.m_ind[i]);

	if (!is_promo(type_num<T>(), pfile.m_type[i]))
		error("wrong type!");
	if (pfile.m_size[i].size() != 0)
		error("wrong dimension!");

	if constexpr (is_Char<T>()) {
		Int temp; fin >> temp; s = (T)temp;
	}
	else if constexpr (is_Int<T>()) {
		fin >> s;
	}
	else if constexpr (is_Doub<T>()) {
		fin >> s;
	}
	else if constexpr (is_Comp<T>()) {
		scanComplex(s, fin);
	}
	else error("unhandled!");
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void load(Vector<T> &v, Str_I varname, Matt_IO pfile)
{
	Long i, n, dim;
	ifstream &fin = pfile.m_in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile.m_ind[i]);

	if (!is_promo(type_num<T>(), pfile.m_type[i]))
		error("wrong type!");
	if (pfile.m_size[i].size() != 1)
		error("wrong dimension!");

	n = pfile.m_size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i) {
		if constexpr (is_Char<T>()) {
			Int temp; fin >> temp;  v[i] = (Char)temp;
		}
		else if constexpr (is_Int<T>() || is_Doub<T>())
			fin >> v[i];
	}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	is_Matrix<Tm>() && (is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void load(Tm &a, Str_I varname, Matt_IO pfile)
{
	Long i, j, dim, m, n;
	ifstream &fin = pfile.m_in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile.m_ind[i]);

	if (!is_promo(type_num<T>(), pfile.m_type[i]))
		error("wrong type!");
	if (pfile.m_size[i].size() != 2)
		error("wrong dimension!");

	m = pfile.m_size[i][0]; n = pfile.m_size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			if constexpr (is_Char<T>()) {
				Int temp; fin >> temp;  a(i, j) = (Char)temp;
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
inline void load(Tm &a, Str_I varname, Matt_IO pfile)
{
	Long i, j, k, dim, m, n, q;
	ifstream &fin = pfile.m_in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile.m_ind[i]);

	if (!is_promo(type_num<T>(), pfile.m_type[i]))
		error("wrong type!");
	if (pfile.m_size[i].size() != 3)
		error("wrong dimension!");
	
	m = pfile.m_size[i][0]; n = pfile.m_size[i][1]; q = pfile.m_size[i][2];
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
