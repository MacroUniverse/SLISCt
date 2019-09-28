// 3d function in spherical harmonic expansion
#pragma once
#include "scalar_arith.h"
#include "cmat.h"

namespace slisc {

template <class T>
class Flm : public Cmat<Vector<T>>
{
private:
    Long m_Nf;
    Long m_Mmin, m_Mmax;
public:
    typedef Cmat<Vector<T>> Base;
    using Base::nrows;
    using Base::ncols;
    using Base::operator();
    Flm();
    Flm(Long_I Lmax, Long_I Nr); // L = 0...NL-1, M = 0
    Flm(Long_I Lmax, Long_I Mmin, Long_I Mmax, Long_I Nr); // L = 0...NL-1, M = Mmin...Mmax
    Long Lmax() const;
    Long Mmin() const;
    Long Mmax() const;
    Long nf() const;
    // check if fr{L,M} exists, for any L, M
    Bool exist(Long_I L, Long_I M) const;
    const Vector<T> &get(Long_I L, Long_I M) const;
    Vector<T> &get(Long_I L, Long_I M);
    const T &operator()(Long_I L, Long_I M, Long_I i) const;
    T &operator()(Long_I L, Long_I M, Long_I i);
    Flm &operator=(const T &rhs);
    template <class T1>
    Flm &operator=(const T1 &rhs);
    ~Flm();
};

template <class T>
Flm<T>::Flm() : m_Nf(0)
{}

template <class T>
Flm<T>::Flm(Long_I Lmax, Long_I Nr) :
    Base(Lmax+1,1), m_Nf(0), m_Mmin(0), m_Mmax(0)
{
    for (Long i = 0; i <= Lmax; ++i) {
        (*this)[i].resize(Nr); ++m_Nf;
    }
}

template <class T>
Flm<T>::Flm(Long_I Lmax, Long_I Mmin, Long_I Mmax, Long_I Nr) :
    Base(Lmax + 1, Mmax - Mmin + 1), m_Nf(0), m_Mmin(Mmin), m_Mmax(Mmax)
{
    for (Long L = 0; L <= Lmax; ++L) {
        for (Long M = max(-L, Mmin); M <= min(L, Mmax); ++M) {
            (*this)(L, M - Mmin).resize(Nr); ++m_Nf;
        }
    }
}

template <class T>
Long Flm<T>::Lmax() const
{
    return nrows() - 1;
}

template <class T>
Long Flm<T>::Mmin() const
{
    return m_Mmin;
}

template <class T>
Long Flm<T>::Mmax() const
{
    return m_Mmax;
}

template <class T>
Long Flm<T>::nf() const
{
    return m_Nf;
}

template <class T>
Bool Flm<T>::exist(Long_I L, Long_I M) const
{
#ifdef SLS_CHECK_BOUNDS
    if (L < 0) SLS_ERR("L < 0 is illegal!");
#endif
    if (L <= Lmax() && M >= m_Mmin && M <= m_Mmax && get(L, M).size() > 0) return true;
    return false;
}

template <class T>
const Vector<T> &Flm<T>::get(Long_I L, Long_I M) const
{
#ifdef SLS_CHECK_BOUNDS
    if (L < 0 || L > Lmax() || M < m_Mmin || M > m_Mmax)
        SLS_ERR("{L,M} out of bound!");
#endif
    return (*this)(L, M - m_Mmin);
}

template <class T>
Vector<T> &Flm<T>::get(Long_I L, Long_I M)
{
#ifdef SLS_CHECK_BOUNDS
    if (!exist(L, M)) SLS_ERR("{L,M} out of bound!");
#endif
    return (*this)(L, M - m_Mmin);
}

template <class T>
const T &Flm<T>::operator()(Long_I L, Long_I M, Long_I i) const
{
    return get(L, M)[i];
}

template <class T>
T &Flm<T>::operator()(Long_I L, Long_I M, Long_I i)
{
    return get(L, M)[i];
}

template <class T>
Flm<T> &Flm<T>::operator=(const T &rhs)
{
    return operator=<T>(rhs);
}

template <class T> template <class T1>
Flm<T> &Flm<T>::operator=(const T1 &rhs)
{
    for (Long L = 0; L <= Lmax(); ++L)
        for (Long M = max(-L, m_Mmin); M <= min(L, m_Mmax); ++M)
            if (exist(L, M))
                get(L, M) = rhs;
    return *this;
}

template <class T>
Flm<T>::~Flm()
{
    for (Long i = 0; i < nrows(); ++i)
        for (Long j = 0; j < ncols(); ++j)
            (*this)(i, j).resize(0);
}


// ==== arithmetics ===

template <class T, class T1, SLS_IF(
    is_scalar<T>() && is_promo<T, T1>()
)>
void operator *=(Flm<T> &flm, const T1 &s)
{
    for (Long L = 0; L <= flm.Lmax(); ++L) {
        for (Long M = flm.Mmin(); M <= flm.Mmax(); ++M) {
            if (flm.exist(L,M))
                flm.get(L,M) *= s;
        }
    }
}

}
