// pure imaginary scalar type and arithmetics

#pragma once
#include "global.h"

namespace slisc {
template <class T>
class ImagNum
{
protected:
	T m_x;
public:
	typedef T value_type;
	ImagNum() {};
	explicit ImagNum(const T &val);
	operator complex<T>();
	T real() const;
	T imag() const;
	void imag(const T &val);
};

template <class T>
inline ImagNum<T>::ImagNum(const T &val) :
	m_x(val)
{}

template <class T>
inline ImagNum<T>::operator complex<T>()
{
	return complex<T>(0, imag());
}

template <class T>
inline T ImagNum<T>::real() const
{
	return 0;
}

template <class T>
inline T ImagNum<T>::imag() const
{
	return m_x;
}

template <class T>
inline void ImagNum<T>::imag(const T &val)
{
	m_x = val;
}

// arithmetic

// real(), imag()

template <class T>
inline T real(const ImagNum<T> &val)
{
	return 0;
}

template <class T>
inline T imag(const ImagNum<T> &val)
{
	return val.imag();
}

// operator+
inline Imag operator+(Imag z1, Imag z2)
{
	return Imag(imag(z1) + imag(z2));
}

inline Comp operator+(Imag z1, Comp z2)
{
	return Comp(real(z2), imag(z1)+imag(z2));
}

inline Comp operator+(Comp z1, Imag z2)
{
	return z2 + z1;
}

// operator*
inline Imag operator*(Imag z, Doub x)
{
	return Imag(z.imag()*x);
}

inline Imag operator*(Doub x, Imag z)
{
	return Imag(z.imag()*x);
}

inline Doub operator*(Imag z1, Imag z2)
{
	return -z1.imag()*z2.imag();
}

inline Comp operator*(Imag z1, Comp z2)
{
	return Comp(-z1.imag()*z2.imag(), z1.imag()*z2.real());
}

inline Comp operator*(Comp z1, Imag z2)
{
	return z2 * z1;
}

} // namespace slisc
