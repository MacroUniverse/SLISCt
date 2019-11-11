// COO sparse matrix
#pragma once
#include "sort.h"
#include "svector.h"
#include "arithmetic.h"

namespace slisc {
template <class T>
class MatCoo : public Vbase<T>
{
private:
    typedef Vbase<T> Base;
    using Base::m_p;
    using Base::m_N;
    Long m_Nr, m_Nc, m_Nnz;
    VecLong m_row, m_col;
    T m_zero = (T)0; // TODO: this could be static inline variable for c++17
    MatCoo() {} // default constructor: uninitialized
public:
    using Base::ptr;
    MatCoo(Long_I Nr, Long_I Nc);
    MatCoo(Long_I Nr, Long_I Nc, Long_I Ncap); // reserve Ncap elements
    MatCoo(const MatCoo &rhs);        // Copy constructor
    Long *row_ptr();
    const Long *row_ptr() const;
    Long *col_ptr();
    const Long *col_ptr() const;
    MatCoo & operator=(const MatCoo &rhs);
    template <class T1, SLS_IF(is_promo<T, T1>())>
    MatCoo & operator=(const MatCoo<T1> &rhs);    // copy assignment (do resize(rhs))
    template <class T1, SLS_IF(is_promo<T, T1>())>
    MatCoo & operator=(const CmatObd<T1> &rhs);
    // inline void operator<<(MatCoo &rhs); // move data and rhs.resize(0, 0); rhs.resize(0)
    void push(const T &s, Long_I i, Long_I j); // add one nonzero element
    void set(const T &s, Long_I i, Long_I j); // change existing element or push new element
    Long n1() const;
    Long n2() const;
    Long size() const; // return m_Nr * m_Nc
    Long nnz() const; // return number of non-zero elements
    Long capacity() const;
    // get single index using double index, return -1 if not found
    Long find(Long_I i, Long_I j) const;
    // reference to an element (element must exist)
    T& ref(Long_I i, Long_I j);
    // double indexing (element need not exist)
    const T& operator()(Long_I i, Long_I j) const;
    T &operator()(Long_I ind); // return element
    const T &operator()(Long_I ind) const;
    Long row(Long_I ind) const; // row index
    Long col(Long_I ind) const; // column index
    void trim(Long_I Nnz); // decrease m_Nnz to Nnz
    void resize(Long_I N); // set m_Nz
    void reserve(Long_I N); // reallocate memory, data will be lost m_Nz = 0
    void reshape(Long_I Nr, Long_I Nc); // change matrix shape
    template <class T1>
    void reserve(const MatCoo<T1> &a);
    template <class T1>
    void reshape(const MatCoo<T1> &a);
    void sort_r(); // sort to row major
};

template <class T>
MatCoo<T>::MatCoo(Long_I Nr, Long_I Nc)
    : m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(0), m_col(0)
{
    m_N = 0;
}

template <class T>
MatCoo<T>::MatCoo(Long_I Nr, Long_I Nc, Long_I Ncap) :
    Base(Ncap), m_Nr(Nr), m_Nc(Nc), m_Nnz(0), m_row(Ncap), m_col(Ncap) {}

template <class T>
MatCoo<T>::MatCoo(const MatCoo<T> &rhs)
{
    SLS_ERR("Copy constructor or move constructor is forbidden, use reference "
         "argument for function input or output, and use \"=\" to copy!");
}

template<class T>
Long * MatCoo<T>::row_ptr()
{
    return m_row.ptr();
}

template <class T>
const Long *MatCoo<T>::row_ptr() const
{
    return m_row.ptr();
}

template<class T>
Long * MatCoo<T>::col_ptr()
{
    return m_col.ptr();
}

template <class T>
const Long *MatCoo<T>::col_ptr() const
{
    return m_col.ptr();
}

template <class T>
MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T> &rhs)
{
    return operator=<T>(rhs);
}

template<class T>
template<class T1, SLS_IF0(is_promo<T, T1>())>
MatCoo<T> & MatCoo<T>::operator=(const CmatObd<T1>& rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
template <class T1, SLS_IF0(is_promo<T, T1>())>
MatCoo<T> & MatCoo<T>::operator=(const MatCoo<T1> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
Long MatCoo<T>::find(Long_I i, Long_I j) const
{
    for (Long n = 0; n < m_Nnz; ++n) {
        if (row(n) == i && col(n) == j)
            return n;
    }
    return -1;
}

template <class T>
T& MatCoo<T>::ref(Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("MatCoo::operator()(i,j): index out of bounds!");
#endif
    Long n = find(i, j);
    if (n < 0)
        SLS_ERR("MatCoo::operator()(i,j): element does not exist!");
    return m_p[n];
}

template <class T>
const T &MatCoo<T>::operator()(Long_I i, Long_I j) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_Nr || j < 0 || j >= m_Nc)
        SLS_ERR("MatCoo::operator()(i,j): index out of bounds!");
#endif
    Long n = find(i, j);
    if (n < 0)
        return m_zero; // never return a (const) reference to a temporary
    return m_p[n];
}

template <class T>
void MatCoo<T>::push(const T &s, Long_I i, Long_I j)
{
#ifdef SLS_CHECK_BOUNDS
    if (i<0 || i>=m_Nr || j<0 || j>=m_Nc)
        SLS_ERR("MatCoo::push(): index out of bounds!");
#endif
#ifdef SLS_CHECK_COO_REPEAT
    Long n;
    for (n = 0; n < m_Nnz; ++n) {
        if (row(n) == i && col(n) == j)
            SLS_ERR("MatCoo::push(s,i,j): element already exists!");
    }
#endif
    if (m_Nnz == m_N) SLS_ERR("MatCoo::add(): out of memory, please reserve!");
    m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
    ++m_Nnz;
}

template <class T>
void MatCoo<T>::set(const T &s, Long_I i, Long_I j)
{
    Long n;
    // change
    for (n = 0; n < m_Nnz; ++n) {
        if (row(n) == i && col(n) == j) {
            m_p[n] = s; m_row[n] = i; m_col[n] = j; return;
        }
    }
    // push
    if (m_Nnz == m_N) SLS_ERR("MatCoo::add(): out of memory, please reserve!");
    m_p[m_Nnz] = s; m_row[m_Nnz] = i; m_col[m_Nnz] = j;
    ++m_Nnz;
}

template <class T>
Long MatCoo<T>::n1() const
{
    return m_Nr;
}

template <class T>
Long MatCoo<T>::n2() const
{
    return m_Nc;
}

template <class T>
Long MatCoo<T>::size() const
{
    return m_Nr * m_Nc;
}

template <class T>
Long MatCoo<T>::nnz() const
{
    return m_Nnz;
}

template <class T>
Long MatCoo<T>::capacity() const
{
    return Base::size();
}

template <class T>
T & MatCoo<T>::operator()(Long_I ind)
{
#ifdef SLS_CHECK_BOUNDS
    if (ind<0 || ind>=m_Nnz)
        SLS_ERR("MatCoo::operator(): subscript out of bounds!");
#endif
    return m_p[ind];
}

template <class T>
const T & MatCoo<T>::operator()(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
    if (ind<0 || ind>=m_Nnz)
        SLS_ERR("MatCoo::operator() const: subscript out of bounds!");
#endif
    return m_p[ind];
}

template <class T>
Long MatCoo<T>::row(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
    if (ind<0 || ind>=m_Nnz)
        SLS_ERR("MatCoo::row() subscript out of bounds");
#endif
    return m_row[ind];
}

template <class T>
Long MatCoo<T>::col(Long_I ind) const
{
#ifdef SLS_CHECK_BOUNDS
    if (ind < 0 || ind >= m_Nnz)
        SLS_ERR("MatCoo::col() subscript out of bounds");
#endif
    return m_col[ind];
}

template <class T>
void MatCoo<T>::trim(Long_I Nnz)
{
#ifdef SLS_CHECK_SHAPE
    if (Nnz < 0)
        SLS_ERR("MatCoo::trim() negative input!");
#endif
    if (Nnz < m_Nnz) m_Nnz = Nnz;
    else if (Nnz > m_Nnz) SLS_ERR("MatCoo::trim(): Nnz > m_Nnz!");
}

template <class T>
void MatCoo<T>::resize(Long_I N)
{
#ifdef SLS_CHECK_BOUNDS
    if (N > m_N)
        SLS_ERR("not enough capacity!");
#endif
    m_Nnz = N;
}

template <class T>
void MatCoo<T>::reserve(Long_I N)
{
    Base::resize(N);
    m_row.resize(N);
    m_col.resize(N);
    m_Nnz = 0;
}

template <class T>
void MatCoo<T>::reshape(Long_I Nr, Long_I Nc)
{
    m_Nr = Nr; m_Nc = Nc;
}

template<class T>
inline void MatCoo<T>::sort_r()
{
    VecLong inds(m_Nnz), order(m_Nnz);
    linspace(order, 0, m_Nnz - 1);
    for (Long i = 0; i < m_Nnz; ++i) {
        inds[i] = m_Nc * m_row[i] + m_col[i];
    }
    sort2(inds, order);
    Svector<T> sli(m_p, m_Nnz);
    reorder(sli, order);
    Svector<Long> sli1;
    sli1.set_size(m_Nnz);
    sli1.set_ptr(m_row.ptr());
    reorder(sli1, order);
    sli1.set_ptr(m_col.ptr());
    reorder(sli1, order);
}

template <class T>
template <class T1>
void MatCoo<T>::reshape(const MatCoo<T1> &a)
{
    reshape(a.n1(), a.n2());
}

template <class T>
template <class T1>
void MatCoo<T>::reserve(const MatCoo<T1> &a)
{
    reserve(a.capacity());
}
} // namespace slisc
