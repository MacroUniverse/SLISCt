// save vectors and matrices defined in "nr3.h" to ".mat" or ".matt" files.
// see README.txt for details
// class types: Doub=0, Comp=1, Int=2, Char=3.

#pragma once
//#define MATFILE_BINARY
//#define MATFILE_DUAL

#include "file.h"

namespace slisc {

// matt file object
class Matt;

// ========== Implementation ============

// Matt class for text mode
class Matt {
public:
    Matt();
    Matt(Str_I fname, Char_I *rw, Int_I precision = 17);
    // delimiter between two numbers, can only be ' ' for now.
    static const Char dlm = ' ';
    Char m_rw; // 'r' for read 'w' for write
    ifstream m_in; // read file
    ofstream m_out; // write file
    Int m_n; // variable numbers
    Str fname; // name of the opened file
    vector<Str> m_name; // variable names
    vector<Int> m_type; // variable types
    vector<vector<Long>> m_size; // variable dimensions
    vector<Long> m_ind; // variable positions (line indices)

    // open a file, return 0 if success
    // return -2 if reading failed (e.g. file is not finished, wrong format)
    Int open(Str_I fname, Char_I *rw, Int_I precision = 17);

	Bool isopen();

    // close a file, if not called, will be called in destructor
    void close();

    // ===== internal functions =====

    // get var names and positions from the end of the file
    // after return, matt.m_ind[i] points to the first matrix element;
    // return 0 if successful, return -1 if failed
    Int get_profile();

    // search a variable by name, return index to m_name[i]
    // return -1 if not found
    Int search(Str_I name);

    // read a scalar from m_in
    template <class T, SLS_IF(is_Char<T>())>
    void read(T &s);

    template <class T, SLS_IF(
        is_Int<T>() || is_Llong<T>() || is_Doub<T>())>
    void read(T &s);

    template <class T, SLS_IF(is_Comp<T>())>
    void read(T &s);

    // write a scalar to m_out
    template <class T, SLS_IF(is_scalar<T>())>
    void write(const T &s);

    ~Matt();
};

// read the next variable after previous delimiter
Long scanInverse(ifstream &fin)
{
    Char c;
    Long ind, i, N;

    ind = fin.tellg();
    for (i = 2; i < 100; ++i) {
        fin.seekg(ind - i); c = fin.get();
        if (c == Matt::dlm)
            break;
    }
    fin >> N;
    fin.seekg(ind - i);
    return N;
}

inline Int Matt::get_profile()
{
    Int i, j, n, temp;
    vector<Long> size;
    Str name;
    ifstream &fin = m_in;

    // read number of variables and their positions
    fin.seekg(0, fin.end);
    Long gmax = fin.tellg();
    m_n = (Int)scanInverse(fin);
    if (m_n < 1)
        return -1;
    m_ind.resize(m_n);
    for (i = 0; i < m_n; ++i) {
        m_ind[i] = scanInverse(fin);
        if (m_ind[i] >= gmax || m_ind[i] < 0)
            return -1;
        if (i > 0 && m_ind[i] <= m_ind[i - 1])
            return -1;
    }

    // loop through each variable
    for (i = 0; i < m_n; ++i) {
        fin.seekg(m_ind[i]);
        // read var name
        fin >> n;
        name.resize(0);
        for (j = 0; j < n; ++j) {
            fin >> temp;
            if (temp <= 0 || temp > 127)
                return -1;
            name.push_back((Char)temp);
        }
        m_name.push_back(name);
        // read var type
        fin >> temp;
        if (temp < 0 || temp > 100)
            return -1;
        m_type.push_back(temp);
        // read var dim
        fin >> n;
        if (n < 0 || n > 10)
            return -1;
        size.resize(0);
        for (j = 0; j < n; ++j) {
            fin >> temp;
            if (temp < 0)
                return -1;
            size.push_back(temp);
        }
        m_size.push_back(size);
        m_ind[i] = fin.tellg();
    }
	return 0;
}

// search variable in file by name
inline Int Matt::search(Str_I name)
{
    for (Int i = 0; i < m_n; ++i)
        if (name == m_name[i])
            return i;
    SLS_WARN("variable name not found: " + name + ", file : " + fname);
    return -1;
}

inline Matt::Matt() {}

inline Matt::Matt(Str_I fname, Char_I * rw, Int_I precision)
{ open(fname, rw, precision); }

Int Matt::open(Str_I fname, Char_I *rw, Int_I precision)
{
	if (isopen())
		close();
    this->fname = fname;
    if (rw[0] == 'w') {

#ifndef SLS_MATT_REPLACE
        if (file_exist(fname)) {
            while (true) {
                if (file_exist(fname)) {
                    SLS_WARN("\n\nfile [" + fname + "] already exist! delete file to continue...\n"
                        "  (define SLS_MATT_REPLACE to replace file by default)\n\n");
                }
                else {
                    break;
                }
                pause(10);
            }
        }
#endif
        m_rw = 'w';
        m_n = 0;
        m_out = ofstream(fname);
        if (!m_out.good())
            SLS_ERR("error: file not created (directory does not exist ?): " + fname);
        m_out.precision(precision);
    }
    else {
        m_rw = 'r';
        m_in = ifstream(fname);
        if (!m_in.good())
            SLS_ERR("error: file not found: " + fname);
        m_in.precision(17);
        return get_profile(); // get var names
    }
    return 0;
}

inline Bool Matt::isopen()
{
	return m_in.is_open() != m_out.is_open();
}

inline void Matt::close()
{
    if (m_rw == 'w') {
        ofstream &fout = m_out;
        // write position of variables
        for (Long i = m_ind.size() - 1; i >= 0; --i)
            fout << m_ind[i] << dlm;
        // write number of variables
        fout << m_n;
        m_out.close();
    }
    else {
        m_in.close();
    }
    m_rw = '\0';
    m_n = 0;
    m_name.clear();
    m_type.clear();
    m_size.clear();
    m_ind.clear();
}

// send one scalar to ofstream
template <class T, SLS_IF0(is_scalar<T>())>
void Matt::write(const T &s)
{
    if (is_real<T>()) {
        m_out << to_num(s) << Matt::dlm;
    }
    else if (is_comp<T>()) {
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

template <class T, SLS_IF0(is_Char<T>())>
void Matt::read(T &s)
{
    Int temp; m_in >> temp; s = (T)temp;
}

template <class T, SLS_IF0(
    is_Int<T>() || is_Llong<T>() || is_Doub<T>())>
void Matt::read(T &s)
{
    m_in >> s;
}

template <class T, SLS_IF0(is_Comp<T>())>
void Matt::read(T &c)
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

inline Matt::~Matt()
{
    if (isopen())
        close();
    else if (m_in.is_open() && m_out.is_open())
        SLS_ERR("unknown!");
}


// save() functions

template <class T, SLS_IF(
    is_Char<T>() || is_Int<T>() || is_Llong<T>() || is_Doub<T>() || is_Comp<T>())>
void save(const T &s, Str_I varname, Matt_IO matt)
{
    Long i, n;
    ofstream &fout = matt.m_out;
    if (!fout.is_open())
        SLS_ERR("matt file not open: " + matt.fname);
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

template <class Tv, class T = contain_type<Tv>, SLS_IF(
    ndims<Tv>() == 1 && is_scalar<T>())>
inline void save(const Tv &v, Str_I varname, Matt_IO matt)
{
    Long i, n;
    ofstream &fout = matt.m_out;
    if (!fout.is_open())
        SLS_ERR("matt file not open!");
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
    (is_Char<T>() || is_Int<T>() || is_Llong<T>() || is_Doub<T>() || is_Comp<T>())
)>
inline void save(const Tm &a, Str_I varname, Matt_IO matt,
    Long_I step1 = 1, Long_I step2 = 1)
{
    Long i, j, m, n;
    ofstream &fout = matt.m_out;
    if (!fout.is_open())
        SLS_ERR("matt file not open!");
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
    m = (a.n1() + step1 - 1) / step1; n = (a.n2() + step2 - 1) / step2;
    fout << 2 << Matt::dlm << m << Matt::dlm << n << Matt::dlm;
    // write matrix data
    for (j = 0; j < n; ++j)
        for (i = 0; i < m; ++i) {
            matt.write(a(step1*i, step2*j));
        }
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense<Tmat>() && ndims<Tmat>() == 3 && is_scalar<T>())>
inline void save(const Tmat &a, Str_I varname, Matt_IO matt,
    Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1)
{
    Long i, j, k, m, n, q;
    ofstream &fout = matt.m_out;
    if (!fout.is_open())
        SLS_ERR("matt file not open!");
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
    m = (a.n1() + step1 - 1) / step1; n = (a.n2() + step2 - 1) / step2;
    q = (a.n3() + step3 - 1) / step3;
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
    if (!fout.is_open())
        SLS_ERR("matt file not open!");
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
        m = (a.n2() + step1 - 1) / step1; n = (a.n3() + step2 - 1) / step2;
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
        m = (a.n3() + step1 - 1) / step1; n = (a.n1() + step2 - 1) / step2;
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
        m = (a.n1() + step1 - 1) / step1; n = (a.n2() + step2 - 1) / step2;
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
    if (!fout.is_open())
        SLS_ERR("matt file not open!");
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
        m = (a.n2() + step1 - 1) / step1; n = (a.n3() + step2 - 1) / step2;
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
        m = (a.n3() + step1 - 1) / step1; n = (a.n1() + step2 - 1) / step2;
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
        m = (a.n1() + step1 - 1) / step1; n = (a.n2() + step2 - 1) / step2;
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

// save one var to one file (replace old file)
template <class T>
inline void save(T &s, Str_I varname, Str_I matt_file)
{
    Matt matt(matt_file, "w");
    save(s, varname, matt);
    matt.close();
}

// for string
void save(Str_I str, Str_I varname, Matt_IO matt)
{
    SvecChar_c sli; sli.set(str.data(), str.size());
    save(sli, varname, matt);
}

// ===== read matt files =====
// return 0 if successful, -1 if variable not found

template <class T, SLS_IF(
    is_Char<T>() || is_Int<T>() || is_Llong<T>() ||
    is_Doub<T>() || is_Comp<T>())>
inline Int load(T &s, Str_I varname, Matt_IO matt)
{
    Int i;
    ifstream &fin = matt.m_in;
    i = matt.search(varname);
    if (i < 0)
        return -1;
    fin.seekg(matt.m_ind[i]);

    if (!is_promo(type_num<T>(), matt.m_type[i]))
        SLS_ERR("wrong type!");
    if (matt.m_size[i].size() != 0)
        SLS_ERR("wrong dimension!");

    matt.read<T>(s);
    return 0;
}

template <class T, SLS_IF(
    is_Char<T>() || is_Int<T>() || is_Llong<T>() ||
    is_Doub<T>() || is_Comp<T>())>
inline Int load(Vector<T> &v, Str_I varname, Matt_IO matt)
{
    Long i, n;
    ifstream &fin = matt.m_in;
    i = matt.search(varname);
    if (i < 0)
        return -1;
    fin.seekg(matt.m_ind[i]);

    if (!is_promo(type_num<T>(), matt.m_type[i]))
        SLS_ERR("wrong type!");
    if (matt.m_size[i].size() != 1)
        SLS_ERR("wrong dimension!");

    n = matt.m_size[i][0]; v.resize(n);
    // read var data
    for (i = 0; i < n; ++i)
        matt.read(v[i]);
    return 0;
}

template <class Tm, class T = contain_type<Tm>, SLS_IF(
    is_dense_mat<Tm>() && (is_Char<T>() || is_Int<T>() || is_Llong<T>() ||
        is_Doub<T>() || is_Comp<T>()))>
inline Int load(Tm &a, Str_I varname, Matt_IO matt)
{
    Long i, j, m, n;
    ifstream &fin = matt.m_in;
    i = matt.search(varname);
    if (i < 0)
        return -1;
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
    return 0;
}

template <class Tmat, class T = contain_type<Tmat>, SLS_IF(
    is_dense<Tmat>() && ndims<Tmat>() == 3 && is_scalar<T>())>
inline Int load(Tmat &a, Str_I varname, Matt_IO matt)
{
    Long i, j, k, m, n, q;
    ifstream &fin = matt.m_in;
    i = matt.search(varname);
    if (i < 0)
        return -1;
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
    return 0;
}

// read one var from one file
template <class T>
inline void load(T &s, Str_I varname, Str_I matt_file)
{
    Matt matt(matt_file, "r");
    load(s, varname, matt);
    matt.close();
}

} // namespace slisc
