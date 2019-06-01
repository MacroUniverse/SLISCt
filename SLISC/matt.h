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
	// delimiter between two numbers, can only be ' ' for now.
	static const Char dlm = ' ';
	char m_rw; // 'r' for read 'w' for write
	ifstream m_in; // read file
	ofstream m_out; // write file
	Int m_n; // variable numbers
	vector<Str> m_name; // variable names
	vector<Int> m_type; // variable types
	vector<vector<Long>> m_size; // variable dimensions
	vector<Long> m_ind; // variable positions (line indices)

	// open a file
	void open(Str fname, Char_I *rw, Int_I precision = 17);

	// close a file
	void close();	

	// ===== internal functions =====

	// get var names and positions from the end of the file
	// after return, matt.m_ind[i] points to the first matrix element;
	void get_profile();

	// search a variable by name, return index to m_name[i]
	Int search(Str_I name);

	// read a complex number from m_in
	void readComplex(Comp &c);

	// read a scalar from m_in
	template <class T, SLS_IF(is_scalar<T>())>
	void read(T &s);

	// write a scalar to m_out
	template <class T, SLS_IF(is_scalar<T>())>
	void write(const T &s);
};

// read the next variable after previous delimiter
Long scanInverse(ifstream &fin)
{
	Char c;
	Long ind, i, N;

	ind = fin.tellg();
	for (i = 2; i < 100; ++i) {
		fin.seekg(ind - i); c = fin.get();
		if (c == Matt::dlm) break;
	}
	fin >> N;
	fin.seekg(ind - i);
	return N;
}

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

// search variable in file by name
inline Int Matt::search(Str_I name)
{
	for (Int i = 0; i < m_n; ++i)
		if (name == m_name[i])
			return i;
	SLS_ERR("variable name not found!");
	return -1;
}

void Matt::open(Str fname, Char_I *rw, Int_I precision)
{
	if (rw[0] == 'w') {
		m_rw = 'w';
		m_n = 0;
		m_out = ofstream(fname);
		m_out.precision(precision);
	}
	else {
		m_rw = 'r';
		m_in = ifstream(fname);
		if (!m_in)
			SLS_ERR("error: file not found: ");
		m_in.precision(17);
		get_profile(); // get var names
	}
}

inline void Matt::close()
{
	Llong i;
	if (m_rw == 'w') {
		ofstream &fout = m_out;
		// write position of variables
		for (i = m_ind.size() - 1; i >= 0; --i)
			fout << m_ind[i] << dlm;
		// write number of variables
		fout << m_n;
		m_out.close();
	}
	else {
		m_in.close();
	}
}

// send one scalar to ofstream
template <class T, SLS_IF0(is_scalar<T>())>
void Matt::write(const T &s)
{
	if constexpr (is_real<T>()) {
		m_out << to_num(s) << Matt::dlm;
	}
	else if constexpr (is_comp<T>()) {
		if (imag(s) == 0)
			m_out << real(s) << Matt::dlm;
		else if (imag(s) < 0)
			m_out << real(s) << imag(s) << "i ";
		else
			m_out << real(s) << '+' << imag(s) << "i ";
	}
	else
		SLS_ERR("unhandled case!");
}

template <class T, SLS_IF0(is_scalar<T>())>
void Matt::read(T &s)
{
	if constexpr (is_Char<T>()) {
		Int temp; m_in >> temp; s = (T)temp;
	}
	else if constexpr (is_Int<T>() || is_Doub<T>())
		m_in >> s;
	else if constexpr (is_Comp<T>())
		readComplex(s);
	else SLS_ERR("unhandled!");
}

inline void Matt::readComplex(Comp &c)
{
	Doub cr = 0, ci = 0;
	Char ch;
	m_in >> cr;
	ch = m_in.get();
	if (ch == Matt::dlm) {
		c = cr; return;
	}
	m_in >> ci;
	if (ch == '-')
		ci *= -1.;
	c = Comp(cr, ci);
	m_in.ignore(100, Matt::dlm);
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
void save(const T &s, Str_I varname, Matt_IO matt)
{
	Long i, n;
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << Matt::dlm;
	}
	// write data type info
	fout << type_num<T>() << Matt::dlm;
	// write dimension info
	fout << 0 << Matt::dlm;
	// write matrix data
	matt.write(s);
}

template <class T, SLS_IF(
	is_Int<T>() || is_Char<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void save(const Vector<T> &v, Str_I varname, Matt_IO matt)
{
	Long i, n;
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << Matt::dlm;
	}
	// write data type info
	fout << type_num<T>() << Matt::dlm;
	// write dimension info
	n = v.size();
	fout << 1 << Matt::dlm << n << Matt::dlm;
	// write matrix data
	for (i = 0; i < n; ++i) {
		matt.write(v[i]);
	}
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	(is_Matrix<Tm>() || is_Cmat<Tm>()) &&
	(is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void save(const Tm &a, Str_I varname, Matt_IO matt,
	Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, m, n;
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << Matt::dlm;
	}
	// write data type info
	fout << type_num<T>() << Matt::dlm;
	// write dimension info
	m = (a.n1() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << Matt::dlm << m << Matt::dlm << n << Matt::dlm;
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			matt.write(a(step1*i, step2*j));
		}
}

template <class T, SLS_IF(
	is_Doub<T>() || is_Comp<T>()
)>
inline void save(const Mat3d<T> &a, Str_I varname, Matt_IO matt,
	Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1)
{
	Long i, j, k, m, n, q;
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << to_num(varname.at(i)) << Matt::dlm;
	}
	// write data type info
	fout << type_num<T>() << Matt::dlm;
	// write dimension info
	m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
	q = (a.dim3() + step3 - 1) / step3;
	fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << q << Matt::dlm;
	// write matrix data
	for (k = 0; k < q; ++k)
	for (j = 0; j < n; ++j)
	for (i = 0; i < m; ++i) {
		matt.write(a(step1*i, step2*j, step3*k));
	}
}

inline void save(Mat3Doub_I &a, Str_I varname, Matt_IO matt,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << Matt::dlm;
	}
	// write data type info
	fout << type_num<Doub>() << Matt::dlm;
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
					fout << a(ind, step1*j, step2*k) << Matt::dlm;
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k)
					fout << a(step2*i, ind, step1*k) << Matt::dlm;
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i)
					fout << a(step1*i, step2*j, ind) << Matt::dlm;
		}
	}
	else
		SLS_ERR("illegal value of xyz");
}

inline void save(Mat3Comp_I &a, Str_I varname, Matt_IO matt,
	Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1)
{
	Long i, j, k, m, n, ind, Nslice{ slice.size() };
	ofstream &fout = matt.m_out;
	++matt.m_n; matt.m_ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << Matt::dlm;
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << Matt::dlm;
	}
	// write data type info
	fout << type_num<Comp>() << Matt::dlm;
	if (xyz == 'x') {
		// write dimension info
		m = (a.dim2() + step1 - 1) / step1; n = (a.dim3() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (i = 0; i < Nslice; ++i) {
			ind = slice[i];
			for (k = 0; k < n; ++k)
				for (j = 0; j < m; ++j)
					matt.write(a(ind, step1*j, step2*k));					
		}
	}
	else if (xyz == 'y') {
		// write dimension info
		m = (a.dim3() + step1 - 1) / step1; n = (a.dim1() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (j = 0; j < Nslice; ++j) {
			ind = slice[j];
			for (i = 0; i < n; ++i)
				for (k = 0; k < m; ++k) {
					matt.write(a(step2*i, ind, step1*k));
				}
		}
	}
	else if (xyz == 'z') {
		// write dimension info
		m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
		fout << 3 << Matt::dlm << m << Matt::dlm << n << Matt::dlm << Nslice << Matt::dlm;
		// write matrix data
		for (k = 0; k < Nslice; ++k) {
			ind = slice[k];
			for (j = 0; j < n; ++j)
				for (i = 0; i < m; ++i)
					matt.write(a(step1*i, step2*j, ind));
		}
	}
	else
		SLS_ERR("illegal value of xyz");
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void load(T &s, Str_I varname, Matt_IO matt)
{
	Int i;
	ifstream &fin = matt.m_in;
	i = matt.search(varname);
	fin.seekg(matt.m_ind[i]);

	if (!is_promo(type_num<T>(), matt.m_type[i]))
		SLS_ERR("wrong type!");
	if (matt.m_size[i].size() != 0)
		SLS_ERR("wrong dimension!");

	matt.read<T>(s);
}

template <class T, SLS_IF(
	is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>()
)>
inline void load(Vector<T> &v, Str_I varname, Matt_IO matt)
{
	Long i, n, dim;
	ifstream &fin = matt.m_in;
	i = matt.search(varname);
	fin.seekg(matt.m_ind[i]);

	if (!is_promo(type_num<T>(), matt.m_type[i]))
		SLS_ERR("wrong type!");
	if (matt.m_size[i].size() != 1)
		SLS_ERR("wrong dimension!");

	n = matt.m_size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i)
		matt.read(v[i]);
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	is_Matrix<Tm>() && (is_Char<T>() || is_Int<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void load(Tm &a, Str_I varname, Matt_IO matt)
{
	Long i, j, dim, m, n;
	ifstream &fin = matt.m_in;
	i = matt.search(varname);
	fin.seekg(matt.m_ind[i]);

	if (!is_promo(type_num<T>(), matt.m_type[i]))
		SLS_ERR("wrong type!");
	if (matt.m_size[i].size() != 2)
		SLS_ERR("wrong dimension!");

	m = matt.m_size[i][0]; n = matt.m_size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i)
			matt.read(a(i, j));
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
	is_Mat3d<Tm>() && (is_Doub<T>() || is_Comp<T>())
)>
inline void load(Tm &a, Str_I varname, Matt_IO matt)
{
	Long i, j, k, dim, m, n, q;
	ifstream &fin = matt.m_in;
	i = matt.search(varname);
	fin.seekg(matt.m_ind[i]);

	if (!is_promo(type_num<T>(), matt.m_type[i]))
		SLS_ERR("wrong type!");
	if (matt.m_size[i].size() != 3)
		SLS_ERR("wrong dimension!");
	
	m = matt.m_size[i][0]; n = matt.m_size[i][1]; q = matt.m_size[i][2];
	a.resize(m, n, q);
	// read var data
	for (k = 0; k < q; ++k)
		for (j = 0; j < n; ++j)
			for (i = 0; i < m; ++i)
				matt.read(a(i, j, k));
}

} // namespace slisc
