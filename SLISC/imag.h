// pure imaginary scalar type and arithmetics

#pragma once
#include "meta.h"

namespace slisc {
template <class T>
class ImagNum
{
protected:
    T m_s;
public:
    typedef T value_type;
    ImagNum() {};
    template <class Ts, SLS_IF(is_promo<T, Ts>())>
    constexpr explicit ImagNum(const Ts &val);
    operator complex<T>() const;
    T real() const;
    T imag() const;
    void imag(const T &val);
};

template <class T>
template <class Ts, SLS_IF0(is_promo<T, Ts>())>
constexpr ImagNum<T>::ImagNum(const Ts &val) :
    m_s(val)
{}

template <class T>
inline ImagNum<T>::operator complex<T>() const
{
    return complex<T>(T(0), imag());
}

template <class T>
inline T ImagNum<T>::real() const
{
    return 0;
}

template <class T>
inline T ImagNum<T>::imag() const
{
    return m_s;
}

template <class T>
inline void ImagNum<T>::imag(const T &val)
{
    m_s = val;
}

// arithmetic

// real(), imag()

template <class T>
inline T real(const ImagNum<T> &val)
{
    return T(0);
}

template <class T>
inline T imag(const ImagNum<T> &val)
{
    return val.imag();
}

template <class T>
inline auto abs(const ImagNum<T> &val)
{
    return val.imag();
}

// operator+
inline Imag operator+(Imag_I z1, Imag_I z2)
{
    return Imag(imag(z1) + imag(z2));
}

inline Comp operator+(Doub_I x, Imag_I y)
{
	return Comp(x, imag(y));
}

inline Comp operator+(Imag_I x, Doub_I y)
{
	return Comp(y, imag(x));
}

inline Comp operator+(Imag_I z1, Comp_I z2)
{
    return Comp(real(z2), imag(z1)+imag(z2));
}

inline Comp operator+(Comp_I z1, Imag_I z2)
{
    return z2 + z1;
}

// operator-
inline Imag operator-(Imag_I z)
{
    return Imag(-z.imag());
}

inline Imag operator-(Imag_I z1, Imag_I z2)
{
    return Imag(z1.imag() - z2.imag());
}

inline Comp operator-(Doub_I x, Imag_I z)
{
    return Comp(x, -z.imag());
}

inline Comp operator-(Imag_I z, Doub_I x)
{
    return Comp(-x, z.imag());
}

inline Comp operator-(Comp_I z1, Imag_I z2)
{
    return Comp(z1.real(), z1.imag() - z2.imag());
}

inline Comp operator-(Imag_I z1, Comp_I z2)
{
    return Comp(-z2.real(), z1.imag() - z2.imag());
}

// operator*
inline Imag operator*(Imag_I z, Doub_I x)
{
    return Imag(z.imag()*x);
}

inline Imag operator*(Doub_I x, Imag_I z)
{
    return Imag(z.imag()*x);
}

inline Doub operator*(Imag_I z1, Imag_I z2)
{
    return -z1.imag()*z2.imag();
}

inline Comp operator*(Imag_I z1, Comp_I z2)
{
    return Comp(-z1.imag()*z2.imag(), z1.imag()*z2.real());
}

inline Comp operator*(Comp_I z1, Imag_I z2)
{
    return Comp(-z2.imag()*z1.imag(), z2.imag()*z1.real());
}

inline Imag operator/(Imag_I z, Doub_I x)
{
    return Imag(z.imag() / x);
}

inline Imag operator/(Doub_I x, Imag_I z)
{
    return Imag(-x / z.imag());
}

inline std::ostream &operator<<(std::ostream &out, Imag_I num)
{
    out << num.imag() << 'i';
    return out;
}

// imaginary unit
constexpr Imag I(1.);

} // namespace slisc
