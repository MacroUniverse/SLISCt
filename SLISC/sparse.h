#pragma once
#include "slisc.h"
// sparse matrix classes

namespace slisc {

template <class T>
class MatCOO 
{
private:
	Long m_Nr, m_Nc, m_N, m_Nmax;
	Vector<T> a;
	Vector<Long> ia, ja;
public:
	MatCOO(): m_Nr(0), m_Nc(0), m_N(0), m_Nmax(0) {}
	MatCOO(Long_I Nr, Long_I Nc, Long_I Nmax): m_Nr(Nr), m_Nc(Nc), m_N(0), m_Nmax(Nmax) {}
	MatCOO(const MatCOO &rhs);		// Copy constructor
	inline MatCOO & operator=(const MatCOO &rhs);	// copy assignment
	inline MatCOO & operator=(const T &rhs);
	inline void operator<<(MatCOO &rhs); // move data and rhs.resize(0, 0)
	inline T& operator()(Long_I i, Long_I j);	// double indexing
	inline const T& operator()(Long_I i, Long_I j) const;
	inline Long nrows() const;
	inline Long ncols() const;
	inline void resize(Long_I Nr, Long_I Nc) // resize (contents preserved, -1: don't change)
	{ m_Nr = Nr; m_Nc = Nc; }
	inline void resize(Long_I N); // resize (contents preserved)
	inline void alloc(Long_I Nmax) { a.resize(Nmax); ia.resize(Nmax); ja.resize(Nmax); }
	template <class T1>
	inline void resize(const MatCOO<T1> &a);
	~MatCOO() {};
};

template <class T>
MatCOO<T>::resize(Long_I N)
{
	if (N > m_Nmax) error("MatCoo::resize(): N too large, reallocate memory!");
	m_N = N;
}

// matrix vector multiplication
template <class T1, class T2>
void mul(Vector<T> &v, const MatCOO<T1> a, const Vector<T2> &v1)
{

}

} // namespace slisc
