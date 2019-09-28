// compile print.cpp to use in debugger
#pragma once
#include "../SLISC/arithmetic.h"
#include "../SLISC/time.h"

void test_print()
{
    using namespace slisc;
    VecInt v(10);
    linspace(v, 1, 10);
    disp(v);

    CmatInt a(4, 4);
    linspace(a, 1, 16);
    disp(a);

    Cmat3Int a3(4, 4, 2);
    linspace(a3, 1, 32);
    disp(a3);

    Cmat4Int a4(4, 4, 2, 2);
    linspace(a4, 1, 64);
    disp(a4);

    McooInt coo(5, 5, 25);
    coo.push(1., 0, 0);
    coo.push(2., 1, 1);
    coo.push(3., 2, 2);
    disp(coo);

    McoohInt cooh(5, 5, 25);
    cooh.push(1., 0, 0);
    cooh.push(2., 1, 1);
    cooh.push(3., 2, 2);
    disp(cooh);

    cout << "set breakpoint here and try print() in debugger!" << endl;
    pause();
}
