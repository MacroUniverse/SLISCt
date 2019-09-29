#pragma once
#include "../SLISC/eig.h"
#include "../SLISC/fedvr.h"

inline double test_fedvr_fun(const double x, const int n)
{
    using namespace slisc;
    Doub y = (x + 1)*(x - 1);
    Doub dx = 2 / (n - 1.);
    for (Int i = 0; i < n - 2; ++i) {
        y *= x - (-1 + dx*(i + 1));
    }
    return y;
}

// test Gauss integration and second derivative matrix
inline void test_gauss()
{
    using namespace slisc;
    Int Ngs = 10;
    VecDoub x0(Ngs), w0(Ngs);

    // test Gauss-Lobatto integration
    // int_{-1}^1 (x+1)(x-1) dx = -4/3

    Int i;
    GaussLobatto(x0, w0);

    VecDoub y(Ngs);
    for (i = 0; i < Ngs; ++i) {
        y[i] = test_fedvr_fun(x0[i], 2);
    }
    Doub res = dot(y, w0);
    if (abs(res + 4. / 3.) > 1e-14) SLS_ERR("failed!");
}

// test second derivative matrix
// y(x) = x^2 - 1, y"(x) = 2
inline void test_D2_mat()
{
    using namespace slisc;
    Int Nfe = 2, Ngs = 10;
    Long Nx = Nfe * (Ngs - 1) - 1;
    VecDoub bounds(Nfe + 1); linspace(bounds, -1., 1.);

    // second derivative matrix
    McooDoub D2s(Nx, Nx);
    VecDoub x(Nx), w(Nx), u(Nx);
    D2_matrix(D2s, x, w, u, bounds, Ngs);
    VecDoub y(Nx);
    for (Long i = 0; i < Nx; ++i) {
        y[i] = test_fedvr_fun(x[i], 2);
    }
    y /= u;
    VecDoub d2y(D2s.n1()); // second derivative
    mul(d2y, D2s, y);
    d2y *= u;
    d2y -= 2;
    if (max_abs(d2y) > 5e-13) SLS_ERR("failed!");
}

// bound states of infinite square well
inline void test_inf_sqr_well()
{
    using namespace slisc;

    // === params ===
    Doub xmin = 0, xmax = 1; // box boundary
    Int Ngs = 10, Nfe = 5; // grid points per finite element (including boundaries)
    // =============

    Long Nx = Nfe * (Ngs - 1) - 1;
    VecDoub bounds(Nfe + 1); linspace(bounds, xmin, xmax);

    VecDoub x(Nx), w(Nx), u(Nx);
    CmatDoub H(Nx, Nx); McooDoub Hs(Nx, Nx); // dense and sparse Hamiltonian
    D2_matrix(Hs, x, w, u, bounds, Ngs);
    H.resize(Hs.n1(), Hs.n2()); H = Hs;
    H *= -0.5; Hs *= -0.5;

    // solve eigen states
    VecDoub eigVal(H.n1()); // eigen values / bound state energies
    CmatDoub eigVec(0, 0); eigVec.resize(H); // eigen vectors / bound states wave functions
    eig_sym(eigVal, eigVec, H);

    // test energies

    VecDoub eng(8);
    for (Int n = 1; n <= 8; ++n)
        eng[n - 1] = SQR(PI) / 2. * SQR(n);
    minus_equals_vv(eng.ptr(), eigVal.ptr(), 8);
    if (max_abs_v(eng.ptr(), 8) > 1e-5) SLS_ERR("failed!");

    // TODO: test wave function using analytical solution
}

// bound states of simple harmonic oscillator
inline void test_SHO()
{
    using namespace slisc;

    // === params ===
    Doub xmin = -10, xmax = 10; // box boundary
    Int Ngs = 10, Nfe = 10; // grid points per finite element (including boundaries)
    // ==============

    Long Nx = Nfe * (Ngs - 1) - 1;
    VecDoub bounds(Nfe + 1); linspace(bounds, xmin, xmax);

    VecDoub x(Nx), w(Nx), u(Nx);
    CmatDoub H(Nx, Nx); McooDoub Hs(Nx, Nx); // dense and sparse Hamiltonian
    D2_matrix(Hs, x, w, u, bounds, Ngs);

    H.resize(Hs.n1(), Hs.n2()); H = Hs;
    H *= -0.5; Hs *= -0.5;

    // add potential to Hamiltonian
    for (Long i = 0; i < x.size(); ++i) {
        H(i, i) += 0.5*SQR(x[i]);
    }

    // solve eigen states
    VecDoub eigVal(H.n1()); // eigen values / bound state energies
    CmatDoub eigVec(0, 0); eigVec.resize(H); // eigen vectors / bound states wave functions
    eig_sym(eigVal, eigVec, H);

    // test energies
    VecDoub Eng(8); linspace(Eng, 0.5, 7.5);
    minus_equals_vv(Eng.ptr(), eigVal.ptr(), 8);
    if (max_abs_v(Eng.ptr(), 8) > 1e-5) SLS_ERR("failed!");

    // TODO: test wave function using analytical solution.
}

void test_fedvr()
{
    test_gauss();
    test_D2_mat();
    test_SHO();
    test_inf_sqr_well();
}
