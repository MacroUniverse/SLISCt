// complex related functions

#pragma once
#include "meta.h"

namespace slisc {

// operator ==
#ifdef SLS_CPP17
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr Bool operator==(const T1 &z1, const T2 &z2)
{
    if constexpr (is_real<T1>()) { // r == c
        return real(z2) == z1 && imag(z2) == 0;
    }
    else if constexpr (is_real<T2>()) { // c == r
        return real(z1) == z2 && imag(z1) == 0;
    }
    else { // c == c
        return real(z1) == real(z2) && imag(z1) == imag(z2);
    }
    return false; // supress compiler warning
}
#else
constexpr Bool operator==(Fcomp_I x, Comp_I y)
{return real(x) == real(y) && imag(x) == imag(y);}
constexpr Bool operator==(Comp_I x, Int_I y)
{return real(x) == y && imag(x) == 0;}
constexpr Bool operator==(Lcomp_I x, Comp_I y)
{return real(x) == real(y) && imag(x) == imag(y);}
#endif

// operator !=
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr Bool operator!=(const T1 &z1, const T2 &z2)
{
    return !(z1 == z2);
}

// operator+=

template <class T, class Tr, SLS_IF(
    is_comp<T>() && is_real<Tr>() && type_num<T>() - type_num<Tr>() > 20
)>
constexpr void operator+=(T &z, const Tr &x)
{
    z.real() += x;
}

//template <class T, class T1, SLS_IF(
//    is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
//)>
//constexpr void operator+=(T &z, const T1 &z1)
//{
//    z += (T)z1;
//}

// operator-=
template <class T, class Tr, SLS_IF(
    is_comp<T>() && is_real<Tr>() && type_num<T>() - type_num<Tr>() > 20
)>
constexpr void operator-=(T &z, const Tr &x)
{
    z.real() -= x;
}

//template <class T, class T1, SLS_IF(
//    is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
//)>
//constexpr void operator-=(T &z, const T1 &z1)
//{
//    z -= (T)z1;
//}

// operator*=
template <class T, class Tr, SLS_IF(
    is_comp<T>() && is_real<Tr>() && type_num<T>() - type_num<Tr>() > 20
)>
constexpr void operator*=(T &z, const Tr &x)
{
    z *= (rm_comp<T>)x;
}

//template <class T, class T1, SLS_IF(
//    is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
//)>
//constexpr void operator*=(T &z, const T1 &z1)
//{
//    z *= (T)z1;
//}

// operator/=
template <class T, class Tr, SLS_IF(
    is_comp<T>() && is_real<Tr>() && type_num<T>() - type_num<Tr>() > 20
)>
constexpr void operator/=(T &z, const Tr &x)
{
    z /= (rm_comp<T>)x;
}

//template <class T, class T1, SLS_IF(
//    is_comp<T>() && is_comp<T1>() && type_num<T>() > type_num<T1>()
//)>
//constexpr void operator/=(T &z, const T1 &z1)
//{
//    z /= (T)z1;
//}

// operator+-*/ between comp and real

#ifdef SLS_CPP17
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr const auto operator+(const T1 &z1, const T2 &z2)
{
    typedef promo_type<T1, T2> Tc;
    if constexpr (is_real<T1>()) { // r + c
        return Tc(z1 + real(z2), imag(z2));
    }
    else if constexpr (is_real<T2>()) { // c + r
        return Tc(real(z1) + z2, imag(z1));
    }
    else if constexpr (type_num<T1>() > type_num<T2>()) {
        // c (large) + c (small)
        return z1 + (Tc)z2;
    }
    // c (small) + c (large)
    return (Tc)z1 + z2;
}
#else
Comp operator+(Comp_I x, Int_I y) { return x + (Doub)y; }
Comp operator+(Int_I x, Comp_I y) { return (Doub)x + y; }
#endif

#ifdef SLS_CPP17
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr const auto operator-(const T1 &z1, const T2 &z2)
{
    typedef promo_type<T1, T2> Tc;
    if constexpr (is_real<T1>()) { // r - c
        return Tc(z1 - real(z2), -imag(z2));
    }
    else if constexpr (is_real<T2>()) { // c - r
        return Tc(real(z1) - z2, imag(z1));
    }
    else if constexpr (type_num<T1>() > type_num<T2>()) {
        // c (large) - c (small)
        return z1 - (Tc)z2;
    }
    // c (small) - c (large)
    return (Tc)z1 - z2;
}
#else
Comp operator-(Comp_I x, Int_I y) { return x - (Doub)y; }
Comp operator-(Int_I x, Comp_I y) { return (Doub)x - y; }
Comp operator-(Fcomp_I x, Comp_I y) { return (Comp)x - y; }
Comp operator-(Comp_I x, Fcomp_I y) { return x - (Comp)y; }
#endif

#ifdef SLS_CPP17
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr const auto operator*(const T1 &z1, const T2 &z2)
{
    typedef promo_type<T1, T2> Tc;
    if constexpr (is_real<T1>()) { // r * c
        return Tc(z1*real(z2), z1*imag(z2));
    }
    else if constexpr (is_real<T2>()) { // c * r
        return Tc(real(z1)*z2, imag(z1)*z2);
    }
    else if constexpr (type_num<T1>() > type_num<T2>()) {
        // c (large) - c (small)
        return z1 * (Tc)z2;
    }
    // c (small) * c (large)
    return (Tc)z1 * z2;
}
#else
Comp operator*(Comp_I x, Int_I y) { return x * (Doub)y; }
Comp operator*(Int_I x, Comp_I y) { return (Doub)x * y; }
#endif

#ifdef SLS_CPP17
template <class T1, class T2, SLS_IF(
    is_scalar<T1>() && is_scalar<T2>() &&
    (is_comp<T1>() || is_comp<T2>()) &&
    !is_same<rm_comp<T1>, rm_comp<T2>>()
)>
constexpr const auto operator/(const T1 &z1, const T2 &z2)
{
    typedef rm_comp<promo_type<T1, T2>> Tr;
    if constexpr (is_real<T1>()) { // r / c
        if constexpr (is_same<complex<Tr>, T2>()) {
            // r (small) / c (large)
            return (Tr)z1 / z2;
        }
        else { // r (large) / c (small)
            return z1 / (complex<Tr>)z2;
        }
    }
    else if constexpr (is_real<T2>()) { // c / r
        if constexpr (is_same<complex<Tr>, T1>()) {
            // c (large) / r (small)
            return z1 / (Tr)z2;
        }
        else { // c (small) / r (large)
            return (complex<Tr>)z1 / z2;
        }
    }
    else if constexpr (type_num<T1>() > type_num<T2>()) {
        // c (large) / c (small)
        return z1 / (complex<Tr>)z2;
    }
    // c (small) / c (large)
    return (complex<Tr>)z1 / z2;
}
#else
Comp operator/(Comp_I x, Int_I y) { return x / (Doub)y; }
Comp operator/(Int_I x, Comp_I y) { return (Doub)x / y; }
Comp operator/(Fcomp_I x, Comp_I y) { return (Comp)x / y; }
Comp operator/(Comp_I x, Fcomp_I y) { return x / (Comp)y; }
Comp operator/(Doub_I x, Fcomp_I y) { return x / (Comp)y; }
Comp operator/(Fcomp_I x, Doub_I y) { return (Comp)x / y; }
#endif

} // namespace slisc
