// row-major matrix container
#pragma once
#include "vector.h"

namespace slisc {
// Matrix Class

// convert MatCoo to dense matrix
template <class T, class T1, SLS_IF(
    is_dense<T>() && is_scalar<T1>() &&
    is_promo<contain_type<T>, T1>()
)>
inline T &coo2dense(T &lhs, const MatCoo<T1> &rhs)
{
#ifdef SLS_CHECK_SHAPE
    if (!shape_cmp(lhs, rhs))
        SLS_ERR("wrong shape!");
#endif
    lhs = contain_type<T>(0);
    for (Long i = 0; i < rhs.nnz(); ++i) {
        lhs(rhs.row(i), rhs.col(i)) = rhs[i];
    }
    return lhs;
}

// convert MatCooH to dense matrix
template <class T, class T1, SLS_IF(
    is_dense<T>() && is_scalar<T1>() &&
    is_promo<contain_type<T>, T1>()
)>
inline T &cooh2dense(T &lhs, const MatCooH<T1> &rhs)
{
    lhs.resize(rhs.n1(), rhs.n2());
    lhs = contain_type<T>(0);
    for (Long i = 0; i < rhs.nnz(); ++i) {
        Long r = rhs.row(i), c = rhs.col(i);
        if (r == c)
            lhs(r, r) = rhs(i);
        else {
            lhs(r, c) = rhs(i);
            lhs(c, r) = conj(rhs(i));
        }
    }
    return lhs;
}

template <class T>
class Matrix : public Vbase<T>
{
private:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_Nr, m_Nc;
    Matrix(); // default constructor: uninitialized
public:
    using Base::ptr;
    using Base::operator();
    using Base::operator=;
    Matrix(Long_I Nr, Long_I Nc);
    Matrix(Long_I Nr, Long_I Nc, const T &s);    //Initialize to constant
    Matrix(Long_I Nr, Long_I Nc, const T *ptr);    // Initialize to array
    Matrix(const Matrix &rhs);        // Copy constructor
    Matrix & operator=(const Matrix &rhs);
    template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
    Matrix & operator=(const Tmat &rhs);    // copy assignment
    template <class T1>
    Matrix & operator=(const MatCoo<T1> &rhs);
    template <class T1>
    Matrix & operator=(const MatCooH<T1> &rhs);
#ifdef _CUSLISC_
    Matrix & operator=(const Gmatrix<T> &rhs) // copy from GPU vector
    { rhs.get(*this); return *this; }
#endif
    void operator<<(Matrix &rhs); // move data and rhs.resize(0, 0)
    T& operator()(Long_I i, Long_I j); // double indexing
    const T& operator()(Long_I i, Long_I j) const;
    Long n1() const;
    Long n2() const;
    void resize(Long_I Nr, Long_I Nc); // resize (contents not preserved)
    template <class T1>
    void resize(const Matrix<T1> &a);
};

template <class T>
Matrix<T>::Matrix() {}

template <class T>
Matrix<T>::Matrix(Long_I Nr, Long_I Nc) : Base(Nr*Nc), m_Nr(Nr), m_Nc(Nc) {}

template <class T>
Matrix<T>::Matrix(Long_I Nr, Long_I Nc, const T &s) : Matrix(Nr, Nc)
{ *this = s; }

template <class T>
Matrix<T>::Matrix(Long_I Nr, Long_I Nc, const T *ptr) : Matrix(Nr, Nc)
{ memcpy(m_p, ptr, m_N*sizeof(T)); }

template <class T>
Matrix<T>::Matrix(const Matrix<T> &rhs) : Matrix()
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T>
inline Matrix<T> & Matrix<T>::operator=(const Matrix<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T> template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
inline Matrix<T> & Matrix<T>::operator=(const Tmat &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T> template <class T1>
inline Matrix<T> & Matrix<T>::operator=(const MatCoo<T1> &rhs)
{
    return coo2dense(*this, rhs);
}

template <class T> template <class T1>
inline Matrix<T> & Matrix<T>::operator=(const MatCooH<T1> &rhs)
{
    return cooh2dense(*this, rhs);
}

template <class T>
inline void Matrix<T>::operator<<(Matrix<T> &rhs)
{
    m_Nr = rhs.m_Nr; m_Nc = rhs.m_Nc;
    rhs.m_Nr = rhs.m_Nc = 0;
    Base::operator<<(rhs);
}

template <class T>
inline T& Matrix<T>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_Nc*i+j];
}

template <class T>
inline const T & Matrix<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_Nc*i+j];
}

template <class T>
inline Long Matrix<T>::n1() const
{ return m_Nr; }

template <class T>
inline Long Matrix<T>::n2() const
{ return m_Nc; }

template <class T>
inline void Matrix<T>::resize(Long_I Nr, Long_I Nc)
{
    if (Nr != m_Nr || Nc != m_Nc) {
        Base::resize(Nr*Nc);
        m_Nr = Nr; m_Nc = Nc;
    }
}

template <class T>
template <class T1>
inline void Matrix<T>::resize(const Matrix<T1> &a)
{ resize(a.n1(), a.n2()); }
}
