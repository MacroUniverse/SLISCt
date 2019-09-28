#pragma once
#include "../SLISC/mparith.h"
#include "../SLISC/disp.h"

// TODO : test MpDoub2MpDoub10(), splus(), sminus(), stimes, sdivide, plus, minus, times, divide
// TODO : convert radix 10 to radix 256

inline void test_mparith()
{
    using std::cout; using std::endl; //using std::conj;
    using namespace slisc;
    MParith a;

    /*MpDoub x1, x2;
    x1.pow = 2;
    x1.x.resize(3); x1.x(0) = 10, x1.x(1) = 1, x1.x(2) = 100;

    MpDoub10 x;
    MpDoub2MpDoub10(x, x1);*/

    //std::string str = a.mpPI(10);
    //cout << str << endl;
    ////if (str != "3.1415926535897932384626433832795028841971693993")
    //    //SLS_ERR("failed!");

    //// test plus
    //Long i, N = 5;
    //VecUchar w(N+30), u(N), v(N);
    //u(0) = 100; v(0) = 110;
    //u(1) = 0; v(1) = 0;
    //for (i = 2; i < N; ++i) {
    //    u(i) = v(i) = 0;
    //}
    //a.mp2str(str, u);
    //cout << "u = " << str << endl;
    //a.mp2str(str, v);
    //cout << "v = " << str << endl;
    //a.mpadd(w, u, v);
    //a.mplsh(w);
    //a.mp2str(str, w);
    //cout << "w = " << str << endl;

    //// test mul
    //a.mpmul(w, u, v);
    //a.mp2str(str, w);
    //cout << "w = " << str << endl;

    //// test sdv
    //Int temp;
    //a.mpsdv(w, u, 100, temp);
    //a.mp2str(str, w);
    //cout << "w = " << str << endl;
}
