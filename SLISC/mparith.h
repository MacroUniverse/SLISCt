#pragma once
#include "slisc.h"
#include "fft.h"

namespace slisc {

// Multi-precision double type in radix 256
struct Mp256
{
    VecUchar x; // one digit before radix
    Long pow; // power
    Long N; // how many bytes of x used
    Uchar &operator[](Long_I ind) { return x[ind]; } // read & write
    Int operator()(Long_I ind) { return x[ind]; } // read only
    // TODO: use only N digits when doing arithmetics
};

typedef const Mp256 Mp256_I;
typedef Mp256 Mp256_O, Mp256_IO;

// Multi-precision double type in radix 10
struct Mp10 {
    VecUchar x; // one digit before radix
    Long pow; // pow
};

// Multiple precision arithmetic operations done on character strings,
// interpreted as radix 256 numbers with the radix point after the first digit.
// Implementations for the simpler operations are listed here.
struct MParith {

    // Adds the unsigned radix 256 numbers u and v, yielding the unsigned result w. To achieve
    // the full available accuracy, the array w must be longer, by one element, than the shorter of
    // the two arrays u and v.
    // u, v only have one digit before radix point, but w has two. Use "mplsh" to fix.
    void mpadd(VecUchar_O w, VecUchar_I u, VecUchar_I v) {
        Int j,n=u.size(),m=v.size(),p=w.size();
        Int n_min=MIN(n,m),p_min=MIN(n_min,p-1);
        Uint ireg=0;
        for (j=p_min-1;j>=0;j--) {
            ireg=u[j]+v[j]+hibyte(ireg);
            w[j+1]=lobyte(ireg);
        }
        w[0]=hibyte(ireg);
        if (p > p_min+1)
            for (j=p_min+1;j<p;j++) w[j]=0;
    }

    // Subtracts the unsigned radix 256 number v from u yielding the unsigned result w.If the
    // result is negative(wraps around), is is returned as 1; otherwise it is returned as 0. To
    // achieve the full available accuracy, the array w must be as long as the shorter of the two
    // arrays u and v.
    void mpsub(Int &is, VecUchar_O w, VecUchar_I u, VecUchar_I v) {
        Int j,n=u.size(),m=v.size(),p=w.size();
        Int n_min=MIN(n,m),p_min=MIN(n_min,p-1);
        Uint ireg=256;
        for (j=p_min-1;j>=0;j--) {
            ireg=255+u[j]-v[j]+hibyte(ireg);
            w[j]=lobyte(ireg);
        }
        is=hibyte(ireg)-1;
        if (p > p_min)
            for (j=p_min;j<p;j++) w[j]=0;
    }

    // Short addition: The integer iv (in the range 0  iv  255) is added to the least significant
    // radix position of unsigned radix 256 number u, yielding result w.To ensure that the result
    // does not require two digits before the radix point, one may first right - shift the operand u
    // so that the first digit is 0, and keep track of multiples of 256 separately.
    void mpsad(VecUchar_O w, VecUchar_I u, const Int iv) {
        Int j,n=u.size(),p=w.size();
        Uint ireg=256*iv;
        for (j=n-1;j>=0;j--) {
            ireg=u[j]+hibyte(ireg);
            if (j+1 < p) w[j+1]=lobyte(ireg);
        }
        w[0]=hibyte(ireg);
        for (j=n+1;j<p;j++) w[j]=0;
    }

    // Short multiplication: The unsigned radix 256 number u is multiplied by the integer iv (in
    // the range 0  iv  255), yielding result w.To ensure that the result does not require two
    // digits before the radix point, one may first right - shift the operand u so that the first digit
    // is 0, and keep track of multiples of 256 separately.
    void mpsmu(VecUchar_O w, VecUchar_I u, const Int iv) {
        Int j,n=u.size(),p=w.size();
        Uint ireg=0;
        for (j=n-1;j>=0;j--) {
            ireg=u[j]*iv+hibyte(ireg);
            if (j < p-1) w[j+1]=lobyte(ireg);
        }
        w[0]=hibyte(ireg);
        for (j=n+1;j<p;j++) w[j]=0;
    }

    // Short division: The unsigned radix 256 number u is divided by the integer iv (in the range
    // 0  iv  255), yielding a quotient w and a remainder ir(with 0  ir  255).To achieve
    // the full available accuracy, the array w must be as long as the array u.
    void mpsdv(VecUchar_O w, VecUchar_I u, const Int iv, Int &ir) {
        Int i,j,n=u.size(),p=w.size(),p_min=MIN(n,p);
        ir=0;
        for (j=0;j<p_min;j++) {
            i=256*ir+u[j];
            w[j]=Uchar(i/iv);
            ir=i % iv;
        }
        if (p > p_min)
            for (j=p_min;j<p;j++) w[j]=0;
    }

    // Ones-complement negate the unsigned radix 256 number u.
    void mpneg(VecUchar_IO u) {
        Int j,n=u.size();
        Uint ireg=256;
        for (j=n-1;j>=0;j--) {
            ireg=255-u[j]+hibyte(ireg);
            u[j]=lobyte(ireg);
        }
    }

    // Move the unsigned radix 256 number v into u.
    // To achieve full accuracy, the array v must be as long as the array u.
    void mpmov(VecUchar_O u, VecUchar_I v) {
        Int j,n=u.size(),m=v.size(),n_min=MIN(n,m);
        for (j=0;j<n_min;j++) u[j]=v[j];
        if (n > n_min)
            for(j=n_min;j<n-1;j++) u[j]=0;
    }

    // Left-shift digits of unsigned radix 256 number u.
    // The final element of the array is set to 0.
    void mplsh(VecUchar_IO u) {
        Int j,n=u.size();
        for (j=0;j<n-1;j++) u[j]=u[j+1];
        u[n-1]=0;
    }

    // get the right-most byte of an Uint
    Uchar lobyte(Uint x) {return (x & 0xff);}
    // get the second right-most byte of an Uint
    Uchar hibyte(Uint x) {return ((x >> 8) & 0xff);}

    // The following, more complicated, methods have discussion and implementation below.
    void mpmul(VecUchar_O w, VecUchar_I u, VecUchar_I v);
    void mpinv(VecUchar_O u, VecUchar_I v);
    void mpdiv(VecUchar_O q, VecUchar_O r, VecUchar_I u, VecUchar_I v);
    void mpsqrt(VecUchar_O w, VecUchar_O u, VecUchar_I v);
    void mp2str(Str_O s, VecUchar_I a, Int_I pow = 1);
    Str mpPI(Int_I np);
};

typedef const Mp10 Mp10_I;
typedef Mp10 Mp10_O, Mp10_IO;

// Mp256 multiplication 
inline void times(Mp256_O x, Mp256_I x1, Mp256_I x2) {
    MParith a;
    a.mpmul(x.x, x1.x, x2.x);
    x.pow = x1.pow + x2.pow + 1;
}

// short plus for Mp256
// add 0-255 to the last digit
inline void splus(Mp256_O x, Mp256_I x1, Int_I n)
{
    MParith a;
    Int N = x1.x.size();
    Int last = a.lobyte(Int(x1.x.end()) + n);
    x.x.resize(N);
    a.mpsad(x.x, x.x, n);
    // no loss
    if (x.x(0) == 0) {
        a.mplsh(x.x); x.pow = x1.pow;
        x.x.end() = last;
    }
    // loss one digit and rounded
    else {
        x.pow = x1.pow + 1;
        if (last > 127)
            x.x.end() += 1; // no need to use splus() again!
    }
}

// short multiplication for Mp256
// has rounding on the last digit
inline void stimes(Mp256_O x, Mp256_I x1, Int_I n) {
    Int N = x1.x.size();
    Int last;
    MParith a;
    x.x.resize(N);
    a.mpsmu(x.x, x1.x, n);
    last = a.lobyte(Int(x1.x.end())*n);
    // no loss
    if (x.x(0) == 0) {
        a.mplsh(x.x); x.pow = x1.pow;
        x.x(N - 1) = last;
    }
    // loss 1 digit, rounded
    else {
        if (last > 127)
            splus(x, x, 1);
        x.pow = x1.pow + 1;
    }
}

// inverse
inline void inv(Mp256_O x, Mp256_I x1) {
    MParith a;
    Mp256 temp; temp.x.resize(x1.x);
    a.mpinv(temp.x, x1.x); x.x = temp.x;
    x.pow = -x1.pow;
}

// Mp256 multiplication 
inline void divide(Mp256_O x, Mp256_I x1, Mp256_I x2) {
    MParith a;
    a.mpmul(x.x, x1.x, x2.x);
    x.pow = x1.pow + x2.pow + 1;
}

// convert 0.256 to Mp256 with any precision
// x.x.resize(N)
inline void Mp2560256(Mp256_O x, Int_I N)
{
    x.pow = -1;
    x.x.resize(N);
    Int i, num = 256, temp;
    for (i = 0; i < N; ++i) {
        num *= 256;
        x.x(i) = temp = num/1000;
        num -= temp*1000;
    }
}

inline void Mp256Inv0256(Mp256_O x, Int_I N)
{
    x.pow = 0;
    x.x.resize(N);
    x.x(0) = 3;    x.x(1) = 232;
    for (Long i = 2; i < N; ++i)
        x.x(i) = 0;
}

// short division for Mp256, 0 < n <= 255
void sdivide(Mp256_O x, Mp256_I x1, Int_I n) {
    Int N = x1.x.size();
    Int rem;
    MParith a;
    x.x.resize(N);
    a.mpsdv(x.x, x1.x, n, rem);
    if (x.x(0) == 0) {
        a.mplsh(x.x); x.pow = x1.pow - 1;
        rem *= 256;
        x.x(N - 1) = rem / n; rem %= n;
    }
    else
        x.pow = x1.pow;

    // rounding for last digit
    if (rem*256/n > 127)
        splus(x, x, 1);
}

// convert radix 10 number to Mp256
inline void createMp256(Mp256_O x, Long_I d, Long_I pow)
{
    // TODO
}

// convert Mp256 to radix 10 number
// TODO : specify precision
inline void Mp2562Mp10(Mp10_O x, Mp256_I x1) {
    Long i;
    Int N = x1.x.size();
    MParith a;
    Mp256 x11; // copy of x1, with one more digit
    x.pow = 0;
    x11.x.resize(N+1);
    for (i = 0; i < N; ++i)
        x11.x(i) = x1.x(i);
    x11.x.end() = 0;
    x11.pow = x1.pow;

    if (x1.pow > 0) {
        while (x11.pow > 0 | x11.x(0) > 9) {
            sdivide(x11, x11, 10); x.pow += 1;
        }
    }
    else if (x1.pow < 0) {
        while (x11.pow < 0) {
            stimes(x11, x11, 10); x.pow -= 1;
        }
    }

    Str str;
    a.mp2str(str, x11.x);
    x.x.resize(str.size()-4);
    x.x(0) = str.at(0) - 48;
    for (i = 2; i < str.size()-3; ++i)
        x.x(i-1) = str.at(i) - 48;
}

// convert from Mp10 to Mp256
// N is precision
inline void Mp102Mp256(Mp256_O x, Mp10_I x1, Long_I N)
{
    x.x.resize(N);
    x.x(0) = x1.x(0);
    // TODO
}


// Uses fast Fourier transform to multiply the unsigned radix 256 integers u[0..n-1] and v[0..m-1],
// yielding a product w[0..n + m - 1].
// u, v have 1 digit before radix point, w has 2
void MParith::mpmul(VecUchar_O w, VecUchar_I u, VecUchar_I v) {
    const Doub RX=256.0;
    Int j,nn=1,n=u.size(),m=v.size(),p=w.size(),n_max=MAX(m,n);
    Doub cy,t;
    while (nn < n_max) nn <<= 1;
    nn <<= 1;
    VecDoub a(nn,0.0),b(nn,0.0);
    for (j=0;j<n;j++) a[j]=u[j];
    for (j=0;j<m;j++) b[j]=v[j];
    realft(a,1);
    realft(b,1);
    b[0] *= a[0];
    b[1] *= a[1];
    for (j=2;j<nn;j+=2) {
        b[j]=(t=b[j])*a[j]-b[j+1]*a[j+1];
        b[j+1]=t*a[j+1]+b[j+1]*a[j];
    }
    realft(b,-1);
    cy=0.0;
    for (j=nn-1;j>=0;j--) {
        t=b[j]/(nn >> 1)+cy+0.5;
        cy=Uint(t/RX);
        b[j]=t-cy*RX;
    }
    if (cy >= RX) throw("cannot happen in mpmul");
    for (j=0;j<p;j++) w[j]=0;
    w[0]=Uchar(cy);
    for (j=1;j<MIN(n+m,p);j++) w[j]=Uchar(b[j-1]);
}

// Character string v[0..m-1] is interpreted as a radix 256 number with the radix point after
// (nonzero)v[0]; u[0..n - 1] is set to the most significant digits of its reciprocal, with the radix
// point after u[0].
void MParith::mpinv(VecUchar_O u, VecUchar_I v) {
    const Int MF=4;
    const Doub BI=1.0/256.0;
    Int i,j,n=u.size(),m=v.size(),mm=MIN(MF,m);
    Doub fu,fv=Doub(v[mm-1]);
    VecUchar s(n+m),r(2*n+m);
    for (j=mm-2;j>=0;j--) {
        fv *= BI;
        fv += v[j];
    }
    fu=1.0/fv;
    for (j=0;j<n;j++) {
        i=Int(fu);
        u[j]=Uchar(i);
        fu=256.0*(fu-i);
    }
    for (;;) {
        mpmul(s,u,v);
        mplsh(s);
        mpneg(s);
        s[0] += Uchar(2);
        mpmul(r,s,u);
        mplsh(r);
        mpmov(u,r);
        for (j=1;j<n-1;j++)
            if (s[j] != 0) break;
        if (j==n-1) return;
    }
}

// Divides unsigned radix 256 integers u[0..n-1] by v[0..m-1] (with m  n required), yielding a
// quotient q[0..n - m] and a remainder r[0..m - 1].
void MParith::mpdiv(VecUchar_O q, VecUchar_O r, VecUchar_I u, VecUchar_I v) {
    const Int MACC=1;
    Int i,is,mm,n=u.size(),m=v.size(),p=r.size(),n_min=MIN(m,p);
    if (m > n) throw("Divisor longer than dividend in mpdiv");
    mm=m+MACC;
    VecUchar s(mm),rr(mm),ss(mm+1),qq(n-m+1),t(n);
    mpinv(s,v);
    mpmul(rr,s,u);
    mpsad(ss,rr,1);
    mplsh(ss);
    mplsh(ss);
    mpmov(qq,ss);
    mpmov(q,qq);
    mpmul(t,qq,v);
    mplsh(t);
    mpsub(is,t,u,t);
    if (is != 0) throw("MACC too small in mpdiv");
    for (i=0;i<n_min;i++) r[i]=t[i+n-m];
    if (p>m) for (i=m;i<p;i++) r[i]=0;
}

// Character string v[0..m-1] is interpreted as a radix 256 number with the radix point after
// v[0]; w[0..n - 1] is set to its square root(radix point after w[0]), and u[0..n - 1] is set to the
// reciprocal thereof(radix point before u[0]).w and u need not be distinct, in which case they
// are set to the square root.
void MParith::mpsqrt(VecUchar_O w, VecUchar_O u, VecUchar_I v) {
    const Int MF=3;
    const Doub BI=1.0/256.0;
    Int i,ir,j,n=u.size(),m=v.size(),mm=MIN(m,MF);
    VecUchar r(2*n),x(n+m),s(2*n+m),t(3*n+m);
    Doub fu,fv=Doub(v[mm-1]);
    for (j=mm-2;j>=0;j--) {
        fv *= BI;
        fv += v[j];
    }
    fu=1.0/sqrt(fv);
    for (j=0;j<n;j++) {
        i=Int(fu);
        u[j]=Uchar(i);
        fu=256.0*(fu-i);
    }
    for (;;) {
        mpmul(r,u,u);
        mplsh(r);
        mpmul(s,r,v);
        mplsh(s);
        mpneg(s);
        s[0] += Uchar(3);
        mpsdv(s,s,2,ir);
        for (j=1;j<n-1;j++) {
            if (s[j] != 0) {
                mpmul(t,s,u);
                mplsh(t);
                mpmov(u,t);
                break;
            }
        }
        if (j<n-1) continue;
        mpmul(x,u,v);
        mplsh(x);
        mpmov(w,x);
        return;
    }
}

// Converts a radix 256 fraction a[0..n-1] (radix point before a[0]) to a decimal fraction represented
// as an ASCII string s[0..m - 1], where m is a returned value.
// NOTE: For simplicity, this routine implements a slow(/ N2) algorithm. Fast
// (/ N lnN), more complicated, radix conversion algorithms do exist.
void MParith::mp2str(Str_O s, VecUchar_I a0, Int_I pow)
{
    const Uint IAZ=48;
    char buffer[4];
    VecUchar a; a = a0;
    Int j,m;

    Int n=a.size();
    m=Int(2.408*n);
    sprintf(buffer,"%d",a[0]);
    s=buffer;
    s += '.';
    mplsh(a);
    for (j=0;j<m;j++) {
        mpsmu(a,a,10);
        s += a[0]+IAZ;
        mplsh(a);
    }
}

// Demonstrate multiple precision routines by calculating and printing the first np bytes of pi.
Str MParith::mpPI(Int_I np) {
    const Uint MACC=2;
    Int ir,j,n=np+MACC;
    Uchar mm;
    Str s;
    VecUchar x(n),y(n),sx(n),sxi(n),z(n),t(n),pi(n),ss(2*n),tt(2*n);
    t[0]=2;
    for (j=1;j<n;j++) t[j]=0;
    mpsqrt(x,x,t);
    mpadd(pi,t,x);
    mplsh(pi);
    mpsqrt(sx,sxi,x);
    mpmov(y,sx);
    for (;;) {
        mpadd(z,sx,sxi);
        mplsh(z);
        mpsdv(x,z,2,ir);
        mpsqrt(sx,sxi,x);
        mpmul(tt,y,sx);
        mplsh(tt);
        mpadd(tt,tt,sxi);
        mplsh(tt);
        x[0]++;
        y[0]++;
        mpinv(ss,y);
        mpmul(y,tt,ss);
        mplsh(y);
        mpmul(tt,x,ss);
        mplsh(tt);
        mpmul(ss,pi,tt);
        mplsh(ss);
        mpmov(pi,ss);
        mm=tt[0]-1;
        for (j=1;j < n-1;j++)
            if (tt[j] != mm) break;
        if (j == n-1) {
            mp2str(s,pi);
            // remove error digits
            s.erase(Int(2.408*np),s.length());
            return s;
        }
    }
}

} // namespace slisc
