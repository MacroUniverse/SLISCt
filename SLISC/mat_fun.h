// functions of square matrix
#pragma once
#include "arithmetic.h"
#include "sparse_arith.h"
#include "eig.h"

namespace slisc {

// out = exp(a*s) for symmetric matrix
void exp_mat_sym(CmatDoub_O out, CmatDoub_I a, Doub_I s)
{
#ifdef SLS_CHECK_SHAPE
    if (a.n1() != a.n2() || !shape_cmp(out, a))
        SLS_ERR("not a square matrix!");
#endif
    Long N = a.n1();
    VecDoub eigVal(N);
    CmatDoub eigVec(N, N);
    eig_sym(eigVal, eigVec, a);
    eigVal *= s;
    exp(eigVal, eigVal);
    CmatDoub temp(N, N);
    mul(temp, eigVec, diag(eigVal));
    trans(eigVec);
    // TODO: using cBLAS multiplication can save the transpose
    mul_gen(out, temp, eigVec);
}

void exp_mat_sym(CmatComp_O out, CmatDoub_I a, Comp_I s)
{
#ifdef SLS_CHECK_SHAPE
    if (a.n1() != a.n2() || !shape_cmp(out, a))
        SLS_ERR("not a square matrix!");
#endif
    Long N = a.n1();
    VecDoub eigVal(N);
    VecComp eigValComp(N);
    CmatDoub eigVec(N, N);
    eig_sym(eigVal, eigVec, a);
    Times(eigValComp, eigVal, s);
    exp(eigValComp, eigValComp);
    CmatComp temp(N, N);
    mul(temp, eigVec, diag(eigValComp));
    trans(eigVec);
    // TODO: using cBLAS multiplication can save the transpose
    mul_gen(out, temp, eigVec);
}

// calculate exp(A*s) * v by U * exp(D*s) * Uh * v;
// where s is a scalar, matrix A = U * D * Uh
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
    ExpA(const CmatDoub &A);
    template <class Tv, class Ts, class Tv1, SLS_IF(
        is_dense_vec<Tv>() && is_scalar<Ts>() && is_dense_vec<Tv1>()
    )>
    void mul(Tv &v, const Ts &s, const Tv1 &v1);
};

inline Long ExpA<Doub>::size() const
{
    return m_diag.size();
}

// A is symmetric matrix, only upper triangle is used
inline ExpA<Doub>::ExpA(const CmatDoub & A)
    : m_diag(A.n1()), m_U(A.n1(), A.n2()),
    m_Uh(A.n1(), A.n2())
{
    Long N = A.n1();
#ifdef SLS_CHECK_SHAPE
    if (N != A.n2())
        SLS_ERR("A must be square matrix!");
#endif
    CmatDoub eigVec(N, N);
    eig_sym(m_diag, eigVec, A);
    m_U = eigVec;
    trans(m_Uh, eigVec);
}

template <class Tv, class Ts, class Tv1, SLS_IF0(
    is_dense_vec<Tv>() && is_scalar<Ts>() && is_dense_vec<Tv1>()
)>
inline void ExpA<Doub>::mul(Tv &v, const Ts &s, const Tv1 &v1)
{
#ifdef SLS_CHECK_SHAPE
    if (size() != v.size())
        SLS_ERR("A must be square matrix!");
#endif
    VecComp u(v1.size());
    slisc::mul(u, m_Uh, v1);
    for (Long i = 0; i < size(); ++i) {
        u[i] *= exp(m_diag[i] * s);
    }
    slisc::mul(v, m_U, u);
}

} // namespace slisc
