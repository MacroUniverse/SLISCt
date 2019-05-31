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

// calculate exp(A*t) * v by U * exp(D*t) * Uh * v;
// where A = U * D * Uh
// A must be symmetric or hermitian
template <class T>
class ExpA
{};

template <>
class ExpA<Doub>
{
public:
	VecDoub m_diag;
	MatDoub m_U;
	MatDoub m_Uh;
	// must be symmetric
	Long size() const;
	ExpA(const CmatDoub &A, Doub_I t);
	template <class Tv, class Tv1, SLS_IF(
		is_dense_vec<Tv>() && is_dense_vec<Tv1>()
	)>
	void mul(Tv &v, const Tv1 &v1);
};

inline Long ExpA<Doub>::size() const
{
	return m_diag.size();
}

// A is symmetric matrix, only upper triangle is used
inline ExpA<Doub>::ExpA(const CmatDoub & A, Doub_I t)
	: m_diag(A.nrows()), m_U(A.nrows(), A.ncols()),
	m_Uh(A.nrows(), A.ncols())
{
	Long N = A.nrows();
#ifdef SLS_CHECK_SHAPE
	if (N != A.ncols())
		SLS_ERR("A must be square matrix!");
#endif
	CmatDoub eigVec(N, N);
	eig_sym(m_diag, eigVec, A);
	m_diag *= t;
	exp(m_diag);
	m_U = eigVec;
	trans(m_Uh, eigVec);
}

template <class Tv, class Tv1, SLS_IF0(
	is_dense_vec<Tv>() && is_dense_vec<Tv1>()
)>
inline void ExpA<Doub>::mul(Tv &v, const Tv1 &v1)
{
#ifdef SLS_CHECK_SHAPE
	if (size() != v.size())
		SLS_ERR("A must be square matrix!");
#endif
	VecComp u(v1.size());
	slisc::mul(u, m_Uh, v1);
	u *= m_diag;
	slisc::mul(v, m_U, u);
}

} // namespace slisc
