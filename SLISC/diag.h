// square diagonal matrix
// mostly a clone a Vector<T>
#pragma once
#include "vector.h"

namespace slisc {
template <class T>
class Diag : public Vector<T>
{
private:
    typedef Vector<T> Base;
    Diag() : Base() {} // default constructor: uninitialized
public:
    using Base::operator();
    using Base::operator=;
    Diag(Long_I N);
    Diag(Long_I N, const T &s);
    Diag(const Vector<T> &v);
    Long size() const;
    Long nnz() const;
    Long n1() const;
    Long n2() const;
    Diag &operator=(const Diag &rhs);
    Diag &operator=(const Vector<T> &rhs);
    T &operator()(Long_I i, Long_I j);
    const T &operator()(Long_I i, Long_I j) const;
};

template <class T>
Diag<T>::Diag(Long_I N) : Base(N) {}

template <class T>
Diag<T>::Diag(Long_I N, const T &s) : Base(N, s) {}

template <class T>
Diag<T>::Diag(const Vector<T> &v) : Base(v.size())
{
    *this = v;
}

template <class T>
Long Diag<T>::size() const
{
    SLS_ERR("use nnz() instead!");
    return 0;
}

template <class T>
Long Diag<T>::nnz() const
{
    return Base::size();
}

template <class T>
Long Diag<T>::n1() const
{
    return Base::size();
}

template <class T>
Long Diag<T>::n2() const
{
    return Base::size();
}

template <class T>
Diag<T> &Diag<T>::operator=(const Diag &rhs)
{
    Base::operator=(rhs); return *this;
}

template <class T>
Diag<T> &Diag<T>::operator=(const Vector<T> &rhs)
{
    Base::operator=(rhs); return *this;
}

template <class T>
T &Diag<T>::operator()(Long_I i, Long_I j)
{
    if (i == j) return (*this)[i];
    return T();
}

template <class T>
const T &Diag<T>::operator()(Long_I i, Long_I j) const
{
    if (i == j) return (*this)[i];
    return T();
}

// convert vector to diagonal matrix
template <class Tv, class T = contain_type<Tv>,
    SLS_IF(is_dense_vec<Tv>())>
const Diag<T> &diag(const Tv &v)
{
    return (Diag<T>&)v;
}
} // namespace slisc
