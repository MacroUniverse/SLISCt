// complex related functions

#pragma once
#include "global.h"

namespace slisc {

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
