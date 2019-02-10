// complex related functions

#pragma once
#include "global.h"
#include "meta.h"

namespace slisc {

// operator ==
template <class T1, class T2, SLS_IF(
	is_comp<T1>() && is_comp<T2>() && !is_same<T1, T2>()
)>
constexpr Bool operator==(const T1 &z1, const T2 &z2)
{
	return real(z1) == real(z2) && imag(z1) == imag(z2);
}

template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && !is_same<rm_comp<T>, Tr>()
)>
constexpr Bool operator==(const T &z, const Tr &x)
{
	return real(z) == x && imag(z) == 0;
}

template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && !is_same<rm_comp<T>, Tr>()
)>
constexpr Bool operator==(const Tr &x, const T &z)
{
	return real(z) == x && imag(z) == 0;
}

// operator!=
template <class T1, class T2, SLS_IF(
	is_comp<T1>() && is_comp<T2>() && !is_same<T1, T2>()
)>
constexpr Bool operator!=(const T1 &z1, const T2 &z2)
{
	return !(z1 == z2);
}

template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && !is_same<rm_comp<T>, Tr>()
)>
constexpr Bool operator!=(const T &z, const Tr &x)
{
	return !(z == x);
}

template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && !is_same<rm_comp<T>, Tr>()
)>
constexpr Bool operator!=(const Tr &x, const T &z)
{
	return !(z == x);
}

// operator+=
template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && is_promo<T, Tr>()
)>
constexpr void operator+=(T &z, const Tr &x)
{
	z.real() += x;
}

template <class T, class T1, SLS_IF(
	is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
)>
constexpr void operator+=(T &z, const T1 &z1)
{
	z += (T)z1;
}

// operator-=
template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && is_promo<T, Tr>()
)>
constexpr void operator-=(T &z, const Tr &x)
{
	z.real() -= x;
}

template <class T, class T1, SLS_IF(
	is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
)>
constexpr void operator-=(T &z, const T1 &z1)
{
	z -= (T)z1;
}

// operator*=
template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && is_promo<T, Tr>()
)>
constexpr void operator*=(T &z, const Tr &x)
{
	z *= (rm_comp<T>)x;
}

template <class T, class T1, SLS_IF(
	is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
)>
constexpr void operator*=(T &z, const T1 &z1)
{
	z *= (T)z1;
}

// operator/=
template <class T, class Tr, SLS_IF(
	is_comp<T>() && is_real<Tr>() && is_promo<T, Tr>()
)>
constexpr void operator/=(T &z, const Tr &x)
{
	z /= (rm_comp<T>)x;
}

template <class T, class T1, SLS_IF(
	is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
)>
constexpr void operator/=(T &z, const T1 &z1)
{
	z /= (T)z1;
}

// TODO: use template for operator+-*/

// for Fcomp
inline const Fcomp operator+(Fcomp_I z, Int_I x) { return z + (Float)x; }
inline const Fcomp operator+(Int_I x, Fcomp_I z) { return z + (Float)x; }
inline const Fcomp operator-(Fcomp_I z, Int_I x) { return z - (Float)x; }
inline const Fcomp operator-(Int_I x, Fcomp_I z) { return (Float)x - z; }
inline const Fcomp operator*(Fcomp_I z, Int_I x) { return z * (Float)x; }
inline const Fcomp operator*(Int_I x, Fcomp_I z) { return (Float)x * z; }
inline const Fcomp operator/(Fcomp_I z, Int_I x) { return z / (Float)x; }
inline const Fcomp operator/(Int_I x, Fcomp_I z) { return (Float)x / z; }

inline const Fcomp operator+(Fcomp_I z, Long_I x) { return z + (Float)x; }
inline const Fcomp operator+(Long_I x, Fcomp_I z) { return z + (Float)x; }
inline const Fcomp operator-(Fcomp_I z, Long_I x) { return z - (Float)x; }
inline const Fcomp operator-(Long_I x, Fcomp_I z) { return (Float)x - z; }
inline const Fcomp operator*(Fcomp_I z, Long_I x) { return z * (Float)x; }
inline const Fcomp operator*(Long_I x, Fcomp_I z) { return (Float)x * z; }
inline const Fcomp operator/(Fcomp_I z, Long_I x) { return z / (Float)x; }
inline const Fcomp operator/(Long_I x, Fcomp_I z) { return (Float)x / z; }

inline const Comp operator+(Fcomp_I z, Doub_I x) { return (Comp)z + x; }
inline const Comp operator+(Doub_I x, Fcomp_I z) { return x + (Comp)z; }
inline const Comp operator-(Fcomp_I z, Doub_I x) { return (Comp)z - x; }
inline const Comp operator-(Doub_I x, Fcomp_I z) { return x - (Comp)z; }
inline const Comp operator*(Fcomp_I z, Doub_I x) { return (Comp)z * x; }
inline const Comp operator*(Doub_I x, Fcomp_I z) { return x * (Comp)z; }
inline const Comp operator/(Fcomp_I z, Doub_I x) { return (Comp)z / x; }
inline const Comp operator/(Doub_I x, Fcomp_I z) { return x / (Comp)z; }

// for Comp
inline const Comp operator+(Comp_I z, Int_I x) { return z + (Doub)x; }
inline const Comp operator+(Int_I x, Comp_I z) { return z + (Doub)x; }
inline const Comp operator-(Comp_I z, Int_I x) { return z - (Doub)x; }
inline const Comp operator-(Int_I x, Comp_I z) { return (Doub)x - z; }
inline const Comp operator*(Comp_I z, Int_I x) { return z * (Doub)x; }
inline const Comp operator*(Int_I x, Comp_I z) { return (Doub)x * z; }
inline const Comp operator/(Comp_I z, Int_I x) { return z / (Doub)x; }
inline const Comp operator/(Int_I x, Comp_I z) { return (Doub)x / z; }

inline const Comp operator+(Comp_I z, Long_I x) { return z + (Doub)x; }
inline const Comp operator+(Long_I x, Comp_I z) { return z + (Doub)x; }
inline const Comp operator-(Comp_I z, Long_I x) { return z - (Doub)x; }
inline const Comp operator-(Long_I x, Comp_I z) { return (Doub)x - z; }
inline const Comp operator*(Comp_I z, Long_I x) { return z * (Doub)x; }
inline const Comp operator*(Long_I x, Comp_I z) { return (Doub)x * z; }
inline const Comp operator/(Comp_I z, Long_I x) { return z / (Doub)x; }
inline const Comp operator/(Long_I x, Comp_I z) { return (Doub)x / z; }

inline const Comp operator+(Comp_I z, Float_I x) { return z + (Doub)x; }
inline const Comp operator+(Float_I x, Comp_I z) { return z + (Doub)x; }
inline const Comp operator-(Comp_I z, Float_I x) { return z - (Doub)x; }
inline const Comp operator-(Float_I x, Comp_I z) { return (Doub)x - z; }
inline const Comp operator*(Comp_I z, Float_I x) { return z * (Doub)x; }
inline const Comp operator*(Float_I x, Comp_I z) { return (Doub)x * z; }
inline const Comp operator/(Comp_I z, Float_I x) { return z / (Doub)x; }
inline const Comp operator/(Float_I x, Comp_I z) { return (Doub)x / z; }

inline const Comp operator+(Comp_I z, Fcomp_I zf) { return z + (Comp)zf; }
inline const Comp operator+(Fcomp_I zf, Comp_I z) { return (Comp)zf + z; }
inline const Comp operator-(Comp_I z, Fcomp_I zf) { return z - (Comp)zf; }
inline const Comp operator-(Fcomp_I zf, Comp_I z) { return (Comp)zf - z; }
inline const Comp operator*(Comp_I z, Fcomp_I zf) { return z * (Comp)zf; }
inline const Comp operator*(Fcomp_I zf, Comp_I z) { return (Comp)zf * z; }
inline const Comp operator/(Comp_I z, Fcomp_I zf) { return z / (Comp)zf; }
inline const Comp operator/(Fcomp_I zf, Comp_I z) { return (Comp)zf / z; }

// for Lcomp
inline const Lcomp operator+(Lcomp_I z, Int_I x) { return z + (Ldoub)x; }
inline const Lcomp operator+(Int_I x, Lcomp_I z) { return z + (Ldoub)x; }
inline const Lcomp operator-(Lcomp_I z, Int_I x) { return z - (Ldoub)x; }
inline const Lcomp operator-(Int_I x, Lcomp_I z) { return (Ldoub)x - z; }
inline const Lcomp operator*(Lcomp_I z, Int_I x) { return z * (Ldoub)x; }
inline const Lcomp operator*(Int_I x, Lcomp_I z) { return (Ldoub)x * z; }
inline const Lcomp operator/(Lcomp_I z, Int_I x) { return z / (Ldoub)x; }
inline const Lcomp operator/(Int_I x, Lcomp_I z) { return (Ldoub)x / z; }

inline const Lcomp operator+(Lcomp_I z, Long_I x) { return z + (Ldoub)x; }
inline const Lcomp operator+(Long_I x, Lcomp_I z) { return z + (Ldoub)x; }
inline const Lcomp operator-(Lcomp_I z, Long_I x) { return z - (Ldoub)x; }
inline const Lcomp operator-(Long_I x, Lcomp_I z) { return (Ldoub)x - z; }
inline const Lcomp operator*(Lcomp_I z, Long_I x) { return z * (Ldoub)x; }
inline const Lcomp operator*(Long_I x, Lcomp_I z) { return (Ldoub)x * z; }
inline const Lcomp operator/(Lcomp_I z, Long_I x) { return z / (Ldoub)x; }
inline const Lcomp operator/(Long_I x, Lcomp_I z) { return (Ldoub)x / z; }

inline const Lcomp operator+(Lcomp_I z, Float_I x) { return z + (Ldoub)x; }
inline const Lcomp operator+(Float_I x, Lcomp_I z) { return (Ldoub)x + z; }
inline const Lcomp operator-(Lcomp_I z, Float_I x) { return z - (Ldoub)x; }
inline const Lcomp operator-(Float_I x, Lcomp_I z) { return (Ldoub)x - z; }
inline const Lcomp operator*(Lcomp_I z, Float_I x) { return z * (Ldoub)x; }
inline const Lcomp operator*(Float_I x, Lcomp_I z) { return(Ldoub)x * z; }
inline const Lcomp operator/(Lcomp_I z, Float_I x) { return z / (Ldoub)x; }
inline const Lcomp operator/(Float_I x, Lcomp_I z) { return (Ldoub)x / z; }

inline const Lcomp operator+(Lcomp_I z, Doub_I x) { return z + (Ldoub)x; }
inline const Lcomp operator+(Doub_I x, Lcomp_I z) { return (Ldoub)x + z; }
inline const Lcomp operator-(Lcomp_I z, Doub_I x) { return z - (Ldoub)x; }
inline const Lcomp operator-(Doub_I x, Lcomp_I z) { return (Ldoub)x - z; }
inline const Lcomp operator*(Lcomp_I z, Doub_I x) { return z * (Ldoub)x; }
inline const Lcomp operator*(Doub_I x, Lcomp_I z) { return (Ldoub)x * z; }
inline const Lcomp operator/(Lcomp_I z, Doub_I x) { return z / (Ldoub)x; }
inline const Lcomp operator/(Doub_I x, Lcomp_I z) { return (Ldoub)x / z; }

inline const Lcomp operator+(Lcomp_I zl, Fcomp_I zf) { return zl + (Lcomp)zf; }
inline const Lcomp operator+(Fcomp_I zf, Lcomp_I zl) { return (Lcomp)zf + zl; }
inline const Lcomp operator-(Lcomp_I zl, Fcomp_I zf) { return zl - (Lcomp)zf; }
inline const Lcomp operator-(Fcomp_I zf, Lcomp_I zl) { return (Lcomp)zf - zl; }
inline const Lcomp operator*(Lcomp_I zl, Fcomp_I zf) { return zl * (Lcomp)zf; }
inline const Lcomp operator*(Fcomp_I zf, Lcomp_I zl) { return (Lcomp)zf * zl; }
inline const Lcomp operator/(Lcomp_I zl, Fcomp_I zf) { return zl / (Lcomp)zf; }
inline const Lcomp operator/(Fcomp_I zf, Lcomp_I zl) { return (Lcomp)zf / zl; }

inline const Lcomp operator+(Lcomp_I zl, Comp_I z) { return zl + (Lcomp)z; }
inline const Lcomp operator+(Comp_I z, Lcomp_I zl) { return (Lcomp)z + zl; }
inline const Lcomp operator-(Lcomp_I zl, Comp_I z) { return zl - (Lcomp)z; }
inline const Lcomp operator-(Comp_I z, Lcomp_I zl) { return (Lcomp)z - zl; }
inline const Lcomp operator*(Lcomp_I zl, Comp_I z) { return zl * (Lcomp)z; }
inline const Lcomp operator*(Comp_I z, Lcomp_I zl) { return (Lcomp)z * zl; }
inline const Lcomp operator/(Lcomp_I zl, Comp_I z) { return zl / (Lcomp)z; }
inline const Lcomp operator/(Comp_I z, Lcomp_I zl) { return (Lcomp)z / zl; }

} // namespace slisc
