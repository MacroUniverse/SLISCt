// fixed-size vector/matrix containers
#pragma once
#include "global.h"

namespace slisc {

// Base Class for vector/matrix
template <class T, Long N>
class FixVec
{
protected:
    T m_data[N]; // pointer to the first element
public:
    typedef T value_type;
    FixVec() {}
    FixVec(const T &s) { *this = s; }
    T* ptr() { return m_data; } // get pointer
    const T* ptr() const { return m_data; }
    static constexpr Long size() { return N; }
    T & operator[](Long_I i);
    const T & operator[](Long_I i) const;
    T & operator()(Long_I i);
    const T & operator()(Long_I i) const;
    T& end(Long_I i = 1);
    const T& end(Long_I i = 1) const;
    FixVec & operator=(const FixVec &rhs);
    template <class T1>
    FixVec & operator=(const FixVec<T1, N> &rhs);
    FixVec & operator=(const T &rhs); // for scalar
    void resize(Long_I n) const; // dummy function, only allow n = N
};

template <class T, Long N>
inline T & FixVec<T, N>::operator[](Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= N)
        SLS_ERR("FixVec subscript out of bounds");
#endif
    return m_data[i];
}

template <class T, Long N>
inline const T & FixVec<T, N>::operator[](Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= N)
        SLS_ERR("FixVec subscript out of bounds");
#endif
    return m_data[i];
}

template <class T, Long N>
inline T & FixVec<T, N>::operator()(Long_I i)
{ return (*this)[i]; }

template <class T, Long N>
inline const T & FixVec<T, N>::operator()(Long_I i) const
{ return (*this)[i]; }

template <class T, Long N>
inline FixVec<T, N> & FixVec<T, N>::operator=(const FixVec<T, N> &rhs)
{ return operator=<T>(rhs); }

template <class T, Long N>
inline FixVec<T, N> & FixVec<T, N>::operator=(const T &rhs)
{
    vecset(m_data, rhs, N);
    return *this;
}

template <class T, Long N>
inline void FixVec<T, N>::resize(Long_I n) const
{
#ifdef SLS_CHECK_SHAPE
    if (n != N) SLS_ERR("cannot resize fixed vector!");
#endif
}

template <class T, Long N> template <class T1>
inline FixVec<T, N> & FixVec<T, N>::operator=(const FixVec<T1, N> &rhs)
{
    veccpy(m_data, rhs.ptr(), N);
    return *this;
}

template <class T, Long N>
inline T & FixVec<T, N>::end(Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i <= 0 || i > N)
        SLS_ERR("index out of bound");
#endif
    return m_data[N-i];
}

template <class T, Long N>
inline const T & FixVec<T, N>::end(Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i <= 0 || i > N)
        SLS_ERR("index out of bound");
#endif
    return m_data[N-i];
}

// Column major Matrix Class

template <class T, Long Nr, Long Nc>
class FixCmat : public FixVec<T, Nr*Nc>
{
private:
    typedef FixVec<T, Nr*Nc> Base;
    using Base::m_data;
public:
    using Base::operator();
    using Base::operator=;
    FixCmat() {}
    FixCmat(const T &s); //Initialize to constant
    FixCmat(const FixCmat &rhs); // Copy constructor
    static constexpr Long n1();
    static constexpr Long n2();
    FixCmat & operator=(const T &rhs);
    FixCmat & operator=(const FixCmat &rhs); // copy assignment
    template <class T1, Long Nr1, Long Nc1>
    FixCmat & operator=(const FixCmat<T1, Nr1, Nc1> &rhs);
    T& operator()(Long_I i, Long_I j);    // double indexing
    const T& operator()(Long_I i, Long_I j) const;
    void resize(Long nr, Long nc) const; // dummy: only for nr = Nr, nc = Nc
};

template <class T, Long Nr, Long Nc>
FixCmat<T, Nr, Nc>::FixCmat(const T &s)
{ *this = s; }

template <class T, Long Nr, Long Nc>
FixCmat<T, Nr, Nc>::FixCmat(const FixCmat<T, Nr, Nc> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference argument for function input or output, and use \"=\" to copy!");
}

template <class T, Long Nr, Long Nc>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const FixCmat<T, Nr, Nc> &rhs)
{
    return operator=<T, Nr, Nc>(rhs);
}

template <class T, Long Nr, Long Nc>
template <class T1, Long Nr1, Long Nc1>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const FixCmat<T1, Nr1, Nc1> &rhs)
{
    Nr = rhs.n1();
    Nc = rhs.n2();
    Base::operator=(rhs);
    return *this;
}

template <class T, Long Nr, Long Nc>
inline FixCmat<T, Nr, Nc> & FixCmat<T, Nr, Nc>::operator=(const T &rhs)
{
    Base::operator=(rhs);
    return *this;
}

template <class T, Long Nr, Long Nc>
inline T & FixCmat<T, Nr, Nc>::operator()(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= Nr || j < 0 || j >= Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_data[i+Nr*j];
}

template <class T, Long Nr, Long Nc>
inline const T & FixCmat<T, Nr, Nc>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= Nr || j < 0 || j >= Nc)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_data[i+Nr*j];
}

template <class T, Long Nr, Long Nc>
constexpr Long FixCmat<T, Nr, Nc>::n1()
{ return Nr; }

template <class T, Long Nr, Long Nc>
constexpr Long FixCmat<T, Nr, Nc>::n2()
{ return Nc; }

template <class T, Long Nr, Long Nc>
inline void FixCmat<T, Nr, Nc>::resize(Long nr, Long nc) const
{
#ifdef SLS_CHECK_SHAPE
    if (nr != Nr || nc != Nc)
        SLS_ERR("cannot resize fixed matrix!");
#endif
}

} // namespace slisc
