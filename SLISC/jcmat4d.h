// matrix slicing with both `step1` and `step2`
// Jcmat4d can be used for both col and row major
#pragma once
#include "vector.h"

namespace slisc {

template <class T>
class Jcmat4d
{
private:
    T *m_p;
    Long m_N;
    Long m_N1, m_N2, m_N3, m_N4;
    Long m_step1, m_step2, m_step3, m_step4;
public:
    Jcmat4d();
    Jcmat4d(const T *ptr, Long_I N1, Long_I N2, Long_I N3, Long_I N4,
        Long_I step1, Long_I step2, Long_I step3, Long_I step4);
    void set(const T *ptr, Long_I N1, Long_I N2, Long_I N3, Long_I N4,
        Long_I step1, Long_I step2, Long_I step3, Long_I step4);

    // === Cmat member functions ===
    Jcmat4d & operator=(const Jcmat4d &rhs);    // copy assignment
    template <class Tmat, SLS_IF(is_dense_mat<Tmat>())>
    Jcmat4d & operator=(const Tmat &rhs);
    Jcmat4d & operator=(const T &rhs);
    //template <class T1>
    //Cmat & operator=(const MatCoo<T1> &rhs);
    //template <class T1>
    //Cmat & operator=(const MatCooH<T1> &rhs);
    T& operator[](Long_I i);    // single indexing (inefficient)
    const T& operator[](Long_I i) const;
    T& operator()(Long_I i);    // same as operator[]
    const T& operator()(Long_I i) const;
    T& operator()(Long_I i, Long_I j, Long_I k, Long_I l);    // double indexing
    const T& operator()(Long_I i, Long_I j, Long_I k, Long_I l) const;
    Long n1() const;
    Long n2() const;
    Long n3() const;
    Long n4() const;
    Long step1() const;
    Long step2() const;
    Long step3() const;
    Long step4() const;
    Long size() const;
};

template <class T>
Jcmat4d<T>::Jcmat4d() {}

template <class T>
Jcmat4d<T>::Jcmat4d(const T *ptr, Long_I N1, Long_I N2, Long_I N3, Long_I N4,
    Long_I step1, Long_I step2, Long_I step3, Long_I step4)
    : m_p((T *)ptr), m_N1(N1), m_N2(N2), m_N3(N3), m_N4(N4), m_N(N1*N2*N3*N4),
    m_step1(step1), m_step2(step2), m_step3(step3), m_step4(step4)
{}

template <class T>
void Jcmat4d<T>::set(const T *ptr, Long_I N1, Long_I N2, Long_I N3, Long_I N4,
    Long_I step1, Long_I step2, Long_I step3, Long_I step4)
{
    m_p = (T *)ptr; m_N1 = N1; m_N2 = N2; m_N3 = N3; m_N4 = N4;
    m_step1 = step1; m_step2 = step2; m_step3 = step3; m_step4 = step4;
    m_N = N1 * N2 * N3 * N4;
}

template <class T>
inline Jcmat4d<T> & Jcmat4d<T>::operator=(const Jcmat4d<T> &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
template <class Tmat, SLS_IF0(is_dense_mat<Tmat>())>
Jcmat4d<T> & Jcmat4d<T>::operator=(const Tmat &rhs)
{
    copy(*this, rhs);
    return *this;
}

template <class T>
inline Jcmat4d<T> & Jcmat4d<T>::operator=(const T &rhs)
{
    T *p = m_p;
    for (Long l = 0; l < m_N4; ++l) {
        for (Long k = 0; k < m_N3; ++k) {
            for (Long j = 0; j < m_N2; ++j) {
                for (Long i = 0; i < m_N1; ++i) {
                    *p = rhs;
                    p += m_step1;
                }
                p += -m_step1 + m_step2;
            }
            p += -m_step2 + m_step3;
        }
        p += -m_step3 + m_step4;
    }    
    return *this;
}

template <class T>
inline T & Jcmat4d<T>::operator[](Long_I i)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N)
        SLS_ERR("out of bound!");
#endif
    Long N1N2 = m_N1 * m_N2;
    Long N1N2N3 = m_N1 * m_N2 * m_N3;
    Long i3 = i % N1N2N3;
    Long i2 = i3 % N1N2;
    return operator()(i2 % m_N1, i2 / m_N1, i3 / N1N2, i / N1N2N3);
}

template <class T>
inline const T & Jcmat4d<T>::operator[](Long_I i) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N)
        SLS_ERR("out of bound!");
#endif
    Long N1N2 = m_N1 * m_N2;
    Long N1N2N3 = m_N1 * m_N2 * m_N3;
    Long i3 = i % N1N2N3;
    Long i2 = i3 % N1N2;
    return operator()(i2 % m_N1, i2 / m_N1, i3 / N1N2, i / N1N2N3);
}

template <class T>
inline T & Jcmat4d<T>::operator()(Long_I i)
{
    return operator[](i);
}

template <class T>
inline const T & Jcmat4d<T>::operator()(Long_I i) const
{
    return operator[](i);
}

template <class T>
inline T & Jcmat4d<T>::operator()(Long_I i, Long_I j, Long_I k, Long_I l)
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 ||
        k < 0 || k > m_N3 || l < 0 || l > m_N4)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * i + m_step2 * j + m_step3 * k + m_step4 * l];
}

template <class T>
inline const T & Jcmat4d<T>::operator()(Long_I i, Long_I j, Long_I k, Long_I l) const
{
#ifdef SLS_CHECK_BOUNDS
    if (i < 0 || i >= m_N1 || j < 0 || j >= m_N2 ||
        k < 0 || k > m_N3 || l < 0 || l > m_N4)
        SLS_ERR("Matrix subscript out of bounds");
#endif
    return m_p[m_step1 * i + m_step2 * j + m_step3 * k + m_step4 * l];
}

template <class T>
inline Long Jcmat4d<T>::n1() const
{
    return m_N1;
}

template <class T>
inline Long Jcmat4d<T>::n2() const
{
    return m_N2;
}

template <class T>
inline Long Jcmat4d<T>::n3() const
{
    return m_N3;
}

template <class T>
inline Long Jcmat4d<T>::n4() const
{
    return m_N4;
}

template <class T>
inline Long Jcmat4d<T>::step1() const
{
    return m_step1;
}

template <class T>
inline Long Jcmat4d<T>::step2() const
{
    return m_step2;
}

template <class T>
inline Long Jcmat4d<T>::step3() const
{
    return m_step3;
}

template <class T>
inline Long Jcmat4d<T>::step4() const
{
    return m_step4;
}

template<class T>
inline Long Jcmat4d<T>::size() const
{
    return m_N;
}
} // namespace slisc
