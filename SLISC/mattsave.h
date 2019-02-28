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

// MATTFile class for text mode
struct MATTFile {
	char rw; // 'r' for read 'w' for write
	std::ifstream in; // read file
	std::ofstream out; // write file
	Int n; // variable numbers
	std::vector<std::string> name; // variable names
	std::vector<Int> type; // variable types
	std::vector<std::vector<Long>> size; // variable dimensions
	std::vector<Long> ind; // variable positions (line indices)
};

MATTFile *mattOpen(std::string fname, Char_I *rw);

void mattClose(MATTFile *pfile);

// mattsave()

void mattsave(Uchar_I s, const std::string &varname, MATTFile *pfile);

void mattsave(Int_I s, const std::string &varname, MATTFile *pfile);

void mattsave(Doub_I s, const std::string &varname, MATTFile *pfile);

void mattsave(Comp_I s, const std::string &varname, MATTFile *pfile);

void mattsave(VecUchar_I &v, const std::string &varname, MATTFile *pfile);

void mattsave(VecInt_I &v, const std::string &varname, MATTFile *pfile);

void mattsave(VecDoub_I &v, const std::string &varname, MATTFile *pfile);

void mattsave(VecComp_I &v, const std::string &varname, MATTFile *pfile);

void mattsave(MatUchar_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1 = 1, Long_I step2 = 1);

void mattsave(MatInt_I &a, const std::string &varname, MATTFile *pfile,
			Long_I step1 = 1, Long_I step2 = 1);

void mattsave(MatDoub_I &a, const std::string &varname, MATTFile *pfile,
			Long_I step1 = 1, Long_I step2 = 1);

void mattsave(MatComp_I &a, const std::string &varname, MATTFile *pfile,
			Long_I step1 = 1, Long_I step2 = 1);

void mattsave(Mat3Doub_I &a, const std::string &varname, MATTFile *pfile,
			Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1);

void mattsave(Mat3Doub_I &a, const std::string &varname, MATTFile *pfile,
			Char_I xyz, const VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1);

void mattsave(Mat3Comp_I &a, const std::string &varname, MATTFile *pfile,
			Long_I step1 = 1, Long_I step2 = 1, Long_I step3 = 1);

void mattsave(Mat3Comp_I &a, const std::string &varname, MATTFile *pfile,
			Char_I xyz, VecInt_I &slice, Long_I step1 = 1, Long_I step2 = 1);

// mattload()

void mattload(Uchar &i, const std::string &varname, MATTFile *pfile);

void mattload(Int &i, const std::string &varname, MATTFile *pfile);

void mattload(Doub &s, const std::string &varname, MATTFile *pfile);

void mattload(Comp &s, const std::string &varname, MATTFile *pfile);

void mattload(VecUchar_O &v, const std::string &varname, MATTFile *pfile);

void mattload(VecInt_O &v, const std::string &varname, MATTFile *pfile);

void mattload(VecDoub_O &v, const std::string &varname, MATTFile *pfile);

void mattload(VecComp_O &v, const std::string &varname, MATTFile *pfile);

void mattload(MatUchar_O &a, const std::string &varname, MATTFile *pfile);

void mattload(MatInt_O &a, const std::string &varname, MATTFile *pfile);

void mattload(MatDoub_O &a, const std::string &varname, MATTFile *pfile);

void mattload(MatComp_O &a, const std::string &varname, MATTFile *pfile);

void mattload(Mat3Doub_O &a, const std::string &varname, MATTFile *pfile);

void mattload(Mat3Comp_O &a, const std::string &varname, MATTFile *pfile);


// ========== Implementation ============

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
	std::string name;
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

MATTFile *mattOpen(std::string fname, Char_I *rw)
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

inline void mattsave(Uchar_I s, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 3 << '\n';
	// write dimension info
	fout << 0 << '\n';
	// write scalar
	fout << (Int)s << '\n';
}

inline void mattsave(Int_I s, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 2 << '\n';
	// write dimension info
	fout << 0 << '\n';
	// write matrix data
	fout << s << '\n';
}

inline void mattsave(Doub_I s, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
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
	// write dimension info
	fout << 0 << '\n';
	// write matrix data
	fout << s << '\n';
}

inline void mattsave(Comp_I s, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
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
	// write dimension info
	fout << 0 << '\n';
	// write matrix data
	if (imag(s) == 0)
		fout << real(s) << '\n';
	else if (imag(s) < 0)
		fout << real(s) << imag(s) << "i\n";
	else
		fout << real(s) << '+' << imag(s) << "i\n";
}

inline void mattsave(VecUchar_I &v, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 3 << '\n';
	// write dimension info
	n = v.size();
	fout << 1 << '\n' << n << '\n';
	// write matrix data
	for (i = 0; i < n; ++i) {
		fout << (Int)v[i] << '\n';
	}
}

inline void mattsave(VecInt_I &v, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 2 << '\n';
	// write dimension info
	n = v.size();
	fout << 1 << '\n' << n << '\n';
	// write matrix data
	for (i = 0; i < n; ++i) {
		fout << v[i] << '\n';
	}
}

inline void mattsave(VecDoub_I &v, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
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
	// write dimension info
	n = v.size();
	fout << 1 << '\n' << n << '\n';
	// write matrix data
	for (i = 0; i < n; ++i) {
		fout << v[i] << '\n';
	}
}

inline void mattsave(VecComp_I &v, const std::string &varname, MATTFile *pfile)
{
	Long i, n;
	Doub cr, ci;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i)
		fout << (Int)varname.at(i) << '\n';
	// write data type info
	fout << 1 << '\n';
	// write dimension info
	n = v.size();
	fout << 1 << '\n' << n << '\n';
	// write matrix data
	for (i = 0; i < n; ++i) {
		cr = real(v[i]); ci = imag(v[i]);
		if (ci == 0)
			fout << cr << '\n';
		else if (ci < 0)
			fout << cr << ci << "i\n";
		else
			fout << cr << '+' << ci << "i\n";
	}
}

inline void mattsave(MatUchar_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2)
{
	Long i, j, m, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 3 << '\n';
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << '\n' << m << '\n' << n << '\n';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			fout << (Int)a(step1*i, step2*j) << '\n';
		}
}

inline void mattsave(MatInt_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2)
{
	Long i, j, m, n;
	std::ofstream &fout = pfile->out;
	++pfile->n; pfile->ind.push_back(fout.tellp());
	// write variable name info
	n = varname.size();
	fout << n << '\n';
	for (i = 0; i < n; ++i) {
		fout << (Int)varname.at(i) << '\n';
	}
	// write data type info
	fout << 2 << '\n';
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << '\n' << m << '\n' << n << '\n';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			fout << a(step1*i, step2*j) << '\n';
		}
}

inline void mattsave(MatDoub_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2)
{
	Long i, j, m, n;
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
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << '\n' << m << '\n' << n << '\n';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			fout << a(step1*i, step2*j) << '\n';
		}
}

inline void mattsave(MatComp_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2)
{
	Long i, j, m, n;
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
	// write dimension info
	m = (a.nrows() + step1 - 1) / step1; n = (a.ncols() + step2 - 1) / step2;
	fout << 2 << '\n' << m << '\n' << n << '\n';
	// write matrix data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			c = a(step1*i, step2*j); cr = real(c); ci = imag(c);
			if (ci == 0)
				fout << cr << '\n';
			else if (ci < 0)
				fout << cr << ci << "i\n";
			else
				fout << cr << '+' << ci << "i\n";
		}
}

inline void mattsave(Mat3Doub_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2, Long_I step3)
{
	Long i, j, k, m, n, q;
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
	// write dimension info
	m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
	q = (a.dim3() + step3 - 1) / step3;
	fout << 3 << '\n' << m << '\n' << n << '\n' << q << '\n';
	// write matrix data
	for (k = 0; k < q; ++k)
	for (j = 0; j < n; ++j)
	for (i = 0; i < m; ++i)
		fout << a(step1*i, step2*j, step3*k) << '\n';
}

inline void mattsave(Mat3Doub_I &a, const std::string &varname, MATTFile *pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1, Long_I step2)
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

inline void mattsave(Mat3Comp_I &a, const std::string &varname, MATTFile *pfile,
	Long_I step1, Long_I step2, Long_I step3)
{
	Long i, j, k, m, n, q;
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
	// write dimension info
	m = (a.dim1() + step1 - 1) / step1; n = (a.dim2() + step2 - 1) / step2;
	q = (a.dim3() + step3 - 1) / step3;
	fout << 3 << '\n' << m << '\n' << n << '\n' << q << '\n';
	// write matrix data
	for (k = 0; k < q; ++k)
		for (j = 0; j < n; ++j)
			for (i = 0; i < m; ++i) {
				c = a(step1*i, step2*j, step3*k); cr = real(c); ci = imag(c);
				if (ci == 0)
					fout << cr << '\n';
				else if (ci < 0)
					fout << cr << ci << "i\n";
				else
					fout << cr << '+' << ci << "i\n";
			}
}

inline void mattsave(Mat3Comp_I &a, const std::string &varname, MATTFile *pfile,
	Char_I xyz, VecInt_I &slice, Long_I step1, Long_I step2)
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
inline Int nameSearch(const std::string &name, MATTFile *pfile)
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

inline void mattload(Uchar &I, const std::string &varname, MATTFile *pfile)
{
	Int i, temp;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	if (pfile->type[i] != 3 || pfile->size[i].size() != 0)
		error("wrong type or dim!");
	// read var data
	fin >> temp; I = Uchar(temp);
}

inline void mattload(Int &I, const std::string &varname, MATTFile *pfile)
{
	Int i;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);
	
	// read var type and dim
	if (pfile->type[i] < 2 || pfile->size[i].size() != 0)
		error("wrong type or dim!");
	// read var data
	fin >> I;
}

inline void mattload(Doub &I, const std::string &varname, MATTFile *pfile)
{
	Int i;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	if (pfile->type[i] == 1 || pfile->size[i].size() != 0)
		error("wrong type or dim!");
	// read var data
	fin >> I;
}

inline void mattload(Comp &I, const std::string &varname, MATTFile *pfile)
{
	Int i;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	if (pfile->size[i].size() != 0)
		error("wrong type or dim!");
	// read var data
	scanComplex(I, fin);
}

inline void mattload(VecUchar_O &v, const std::string &varname, MATTFile *pfile)
{
	Long i, n, dim, temp;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] != 3 || dim != 1)
		error("wrong type or dim!");
	n = pfile->size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i) {
		fin >> temp;  v[i] = (Uchar)temp;
	}
}

inline void mattload(VecInt_O &v, const std::string &varname, MATTFile *pfile)
{
	Long i, dim, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] < 2 || dim != 1)
		error("wrong type or dim!");
	n = pfile->size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i)
		fin >> v[i];
}

inline void mattload(VecDoub_O &v, const std::string &varname, MATTFile *pfile)
{
	Long i, dim, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] == 1 || dim != 1)
		error("wrong type or dim!");
	n = pfile->size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i)
		fin >> v[i];
}

inline void mattload(VecComp_O &v, const std::string &varname, MATTFile *pfile)
{
	Long i, dim, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (dim != 1)
		error("wrong type or dim!");
	n = pfile->size[i][0]; v.resize(n);
	// read var data
	for (i = 0; i < n; ++i)
		scanComplex(v[i], fin);
}

inline void mattload(MatUchar_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, dim, m, n, temp;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] != 3 || dim != 2)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i) {
			fin >> temp;  a(i, j) = (Uchar)temp;
		}
}

inline void mattload(MatInt_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, dim, m, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] < 2 || dim != 2)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i)
			fin >> a(i, j);
}

inline void mattload(MatDoub_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, dim, m, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] == 1 || dim != 2)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i)
			fin >> a(i, j);
}

inline void mattload(MatComp_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, dim, m, n;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (dim != 2)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; a.resize(m, n);
	// read var data
	for (j = 0; j < n; ++j)
		for (i = 0; i < m; ++i)
			scanComplex(a(i, j), fin);
}

inline void mattload(Mat3Doub_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, k, dim, m, n, q;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (pfile->type[i] == 1 || dim != 3)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; q = pfile->size[i][2];
	a.resize(m, n, q);
	// read var data
	for (k = 0; k < q; ++k)
		for (j = 0; j < n; ++j)
			for (i = 0; i < m; ++i)
				fin >> a(i, j, k);
}

inline void mattload(Mat3Comp_O &a, const std::string &varname, MATTFile *pfile)
{
	Long i, j, k, dim, m, n, q;
	std::ifstream &fin = pfile->in;
	i = nameSearch(varname, pfile);
	fin.seekg(pfile->ind[i]);

	// read var type and dim
	dim = pfile->size[i].size();
	if (dim != 3)
		error("wrong type or dim!");
	m = pfile->size[i][0]; n = pfile->size[i][1]; q = pfile->size[i][2];
	a.resize(m, n, q);
	// read var data
	for (k = 0; k < q; ++k)
		for (j = 0; j < n; ++j)
			for (i = 0; i < m; ++i)
				scanComplex(a(i, j, k), fin);
}

} // namespace slisc
