#include "nr3plus.h"

struct LUdcmp
{
    Int n;
    MatDoub lu;
    VecInt indx;
    Doub d;
    LUdcmp(MatDoub_I &a);
    void solve(VecDoub_I &b, VecDoub_O &x);
    void solve(MatDoub_I &b, MatDoub_O &x);
    void inverse(MatDoub_O &ainv);
    Doub det();
    void mprove(VecDoub_I &b, VecDoub_IO &x);
    MatDoub_I &aref;
};
