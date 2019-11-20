// band diagonal matrix class
#include "cmat.h"

namespace slisc {

template <class T, SLS_IF(is_scalar<T>())>
void mat2band(Cmat<T> &b, const Cmat<T> &a, Long_I Nup, Long_I Nlow);

template <class T, SLS_IF(is_scalar<T>())>
void mat2band(Matrix<T> &b, const Matrix<T> &a, Long_I Nup, Long_I Nlow);

template <class T, SLS_IF(is_scalar<T>())>
void band2mat(Cmat<T> &a, const Cmat<T> &b, Long_I Nup, Long_I Nlow);

template <class T, SLS_IF(is_scalar<T>())>
void band2mat(Matrix<T> &a, const Matrix<T> &b, Long_I Nup, Long_I Nlow);

template <class T>
class Band
{
public:
    Long m_N1;
    Long m_N2;
    Long m_Nup;
    Long m_Nlow;
    Cmat<T> m_a;

    Band(Long_I N1, Long_I N2, Long_I Nup, Long_I Nlow);
    Band(const Cmat<T> &a, Long_I Nup, Long_I Nlow);

    T * ptr();
    const T * ptr() const;
    Long n1() const;
    Long n2() const;
    Long nup() const;
    Long nlow() const;
};

template<class T>
inline Band<T>::Band(Long_I N1, Long_I N2, Long_I Nup, Long_I Nlow):
    m_N1(N1), m_N2(N2), m_Nup(Nup), m_Nlow(Nlow), m_a(Nup+Nlow+1, N2)
{}

template<class T>
inline T * Band<T>::ptr()
{
    return m_a.ptr();
}

template<class T>
inline const T * Band<T>::ptr() const
{
    return m_a.ptr();
}

template<class T>
inline Long Band<T>::n1() const
{
    return m_N1;
}

template<class T>
inline Long Band<T>::n2() const
{
    return m_N2;
}

template<class T>
inline Long Band<T>::nup() const
{
    return m_Nup;
}

template<class T>
inline Long Band<T>::nlow() const
{
    return m_Nlow;
}

template<class T>
inline Band<T>::Band(const Cmat<T> &a, Long_I Nup, Long_I Nlow):
    m_N1(a.n1()), m_N2(a.n2()), m_Nup(Nup), m_Nlow(Nlow), m_a(Nup+Nlow+1, a.n2())
{
    mat2band(m_a, a, m_Nup, m_Nlow);
}

} // namespace slisc
