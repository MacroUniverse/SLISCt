// functions of square matrix
#pragma once
#include "arithmetic.h"
#include "sparse_arith.h"
#include "eig.h"

namespace slisc {

// out = exp(a*t)
void exp_mat_sym(CmatDoub_O out, CmatDoub_I a, Doub_I t)
{
#ifdef SLS_CHECK_SHAPE
	if (a.nrows() != a.ncols() || !shape_cmp(out, a))
		SLS_ERR("not a square matrix!");
#endif
	Long N = a.nrows();
	VecDoub eigVal(N); CmatDoub eigVec(N, N);
	eig_sym(eigVal, eigVec, a);
	eigVal *= t;
	exp(eigVal, eigVal);
	CmatDoub temp(N, N);
	mul(temp, eigVec, diag(eigVal));
	trans(eigVec);
	mul(out, temp, eigVec);
}

// calculate exp(A) * v by U * D * UH;
// A must be symmetric or hermitian
// diagonalize a matrix A = U * D * UH
template <class T>
class ExpA
{};

template <>
class ExpA<Doub>
{
public:
	VecDoub m_expD;
	MatDoub m_U;
	MatDoub m_Uh;
	// must be symmetric
	Long size() const;
	ExpA(const CmatDoub &A);
	VecComp mul(VecComp_I v);
};

inline Long ExpA<Doub>::size() const
{
	return m_expD.size();
}

inline ExpA<Doub>::ExpA(const CmatDoub & A)
	: m_expD(A.nrows()), m_U(A.nrows(), A.ncols()),
	m_Uh(A.nrows(), A.ncols())
{
	Long N = A.nrows();
#ifdef SLS_CHECK_SHAPE
	if (N != A.ncols())
		SLS_ERR("A must be square matrix!");
#endif
	CmatDoub eigVec(N, N);
	eig_sym(m_expD, eigVec, A);
	exp(m_expD);
	m_U = eigVec;
	trans(m_Uh, eigVec);
}

inline VecComp ExpA<Doub>::mul(VecComp_I v)
{
#ifdef SLS_CHECK_SHAPE
	if (size() != v.size())
		SLS_ERR("A must be square matrix!");
#endif
	VecComp u(v.size());
	slisc::mul(u, m_Uh, v);
}

} // namespace slisc
