// sparse matrix containers
#pragma once
#include "matcoo.h"

#ifndef NDEBUG
// make sure (i,j) element doesn't exist when using MatCoo<T>::push(s,i,j) 
#define SLS_CHECK_COO_REPEAT
#endif

namespace slisc {

// sparse Hermitian / symmetric
// only stores the upper triangle
// nnz() is the actual # of non-zero elem. stored
template <class T>
class MatCooH : public MatCoo<T>
{
private:
    typedef MatCoo<T> Base;
    MatCooH() : Base() {} // default constructor : uninitialized
public:
    MatCooH(Long_I Nr, Long_I Nc);
    MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz);
    using Base::operator();
    T &ref(Long_I i, Long_I j); // reference to an element
    const T operator()(Long_I i, Long_I j) const; // double indexing (element need not exist)
    void push(const T &s, Long_I i, Long_I j); // add one nonzero element
    void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
    void reshape(Long_I Nr, Long_I Nc); // change matrix shape
    template <class T1>
    void reshape(const MatCoo<T1> &a);
    template <class T1>
    MatCooH &operator=(const MatCooH<T1> &rhs);
};

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc) : Base(Nr, Nc)
{
#ifdef SLS_CHECK_SHAPE
    if (Nr != Nc) SLS_ERR("must be square matrix!");
#endif
}

template <class T>
MatCooH<T>::MatCooH(Long_I Nr, Long_I Nc, Long_I Nnz) : Base(Nr, Nc, Nnz)
{
#ifdef SLS_CHECK_SHAPE
    if (Nr != Nc) SLS_ERR("must be square matrix!");
#endif
}

// cannot return a const reference since conj() might create a temporary
template <class T>
const T MatCooH<T>::operator()(Long_I i, Long_I j) const
{
    if (i > j) {
        return CONJ(Base::operator()(j, i));
    }        
    return Base::operator()(i, j);
}

template <class T>
T &MatCooH<T>::ref(Long_I i, Long_I j)
{
    if (i > j) {
        SLS_ERR("lower triangle is empty!");
        return Base::ref(j, i);
    }
    else
        return Base::ref(i, j);
}

template <class T>
void MatCooH<T>::push(const T &s, Long_I i, Long_I j)
{
    if (i > j)
        Base::push(CONJ(s), j, i);
    else
        Base::push(s, i, j);
}

template <class T>
void MatCooH<T>::set(const T &s, Long_I i, Long_I j)
{
    if (i > j)
        Base::set(s, j, i);
    else
        Base::set(s, i, j);
}

template <class T>
void MatCooH<T>::reshape(Long_I Nr, Long_I Nc)
{
#ifdef SLS_CHECK_SHAPE
    if (Nr != Nc) SLS_ERR("must be a square matrix!");
#endif
    Base::reshape(Nr, Nc);
}

template <class T> template <class T1>
void MatCooH<T>::reshape(const MatCoo<T1> &a)
{
#ifdef SLS_CHECK_SHAPE
    if (a.n1() != a.n2())
        SLS_ERR("a is not square matrix!");
#endif
    reshape(a.n1());
}

template <class T> template <class T1>
MatCooH<T> &MatCooH<T>::operator=(const MatCooH<T1> &rhs)
{
    Base(*this) = Base(rhs);
}

} // namespace slisc
