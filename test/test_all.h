// do all available tests of SLISC
// TODO: test mparith.h not finished
#pragma once

#include "test_meta.h"
void print(slisc::Cmat3Comp_I v, slisc::Long_I i, slisc::Long_I n1, slisc::Long_I j, slisc::Long_I n2, slisc::Long_I k, slisc::Long_I n3);
#include "test_scalar_arith.h"
#include "test_imag.h"
#include "test_dense.h"
#include "test_cmat4d.h"
#include "test_slice.h"
#include "test_arithmetic.h"
#include "test_fixsize.h"
#include "test_sparse.h"
#include "test_cmatobd.h"
#include "test_interp1.h"
#include "test_fft.h"
#include "test_random.h"
#include "test_sort.h"
//#include "test_eigen_basics.h"
//#include "test_eigen_linsolve.h"
//#include "test_eigen_fft.h"
#include "test_time.h"
#ifdef SLS_USE_GSL
#include "test_ylm.h"
#endif
#include "test_coulomb.h"
#include "test_input.h"
#include "test_disp.h"
#include "test_print.h"

#if defined(SLS_USE_MKL) || defined(SLS_USE_LAPACKE) && defined(SLS_USE_CBLAS)
#include "test_lin_eq.h"
#include "test_eig.h"
#include "test_mat_fun.h"
#include "test_expokit.h"
#include "test_fedvr.h"
#endif

#include "test_except.h"
#include "test_mattsave.h"
#include "test_anglib.h"
#ifdef SLS_USE_GSL
#include "test_gsl.h"
#endif
#include "test_unicode.h"
#include "test_omp.h"
#include "test_search.h"
#include "test_tree.h"

// #include "test/test_mparith.h"

inline void test_all()
{
    using slisc::cout; using slisc::endl;
    using slisc::Input;

    cout << "test_meta()" << endl;
    test_meta();
    cout << "test_scalar_arith()" << endl;
    test_scalar_arith();
    cout << "test_dense()" << endl;
    test_dense();
    cout << "test_cmat4d()" << endl;
    test_cmat4d();
    cout << "test_slice()" << endl;
    test_slice();
    cout << "test_arithmetic()" << endl;
    test_arithmetic();
    cout << "test_imag()" << endl;
    test_imag();
    cout << "test_fixsize()" << endl;
    test_fixsize();
    cout << "test_sparse()" << endl;
    test_sparse();
    cout << "test_cmatobd()" << endl;
    test_cmatobd();
#if defined(SLS_USE_MKL) || defined(SLS_USE_LAPACKE) && defined(SLS_USE_CBLAS)
	cout << "test_lin_eq()" << endl;
    test_lin_eq();
    cout << "test_eig()" << endl;
    test_eig();
    cout << "test_mat_fun()" << endl;
    test_mat_fun();
    cout << "test_expokit()" << endl;
    test_expokit();
    cout << "test_fedvr()" << endl;
    test_fedvr();
#endif
    cout << "test_interp1()" << endl;
    test_interp1();
    cout << "test_fft()" << endl;
    test_fft();
    cout << "test_rand()" << endl;
    test_random();
    cout << "test_sort()" << endl;
    test_sort();
    cout << "test_search()" << endl;
    test_search();
    cout << "test_tree()" << endl;
    test_tree();
    cout << "test_time()" << endl;
    test_time();
#ifdef SLS_USE_GSL
    cout << "test_ylm()" << endl;
    test_ylm();
#endif
    cout << "test_coulomb()" << endl;
    test_coulomb();
    cout << "test_mattsave()" << endl;
    test_mattsave();
    cout << "test_anglib()" << endl;
    test_anglib();
#ifdef SLS_USE_GSL
    cout << "test_gsl()" << endl;
    test_gsl();
#endif
    cout << "test_unicode()" << endl;
    test_unicode();
    cout << "test_omp()" << endl;
    test_omp();

    // eigen
    /*cout << "test_eigen_basics()" << endl;
    test_eigen_basics();
    cout << "test_eigen_linsolve()" << endl;
    test_eigen_linsolve();
    cout << "test_eigen_fft()" << endl;
    test_eigen_fft();*/

    cout << "auto test successful!\n" << endl;
    Input inp;

    if (inp.iBool("run optional tests?")) {
        if (inp.iBool("test_disp() ?")) {
            test_disp();
        }
        if (inp.iBool("test_print() ?")) {
            test_print();
        }
        if (inp.iBool("test_except() ?")) {
            test_except();
        }
    }

    //cout << "test_mparith()" << endl;
    //test_mparith();
    cout << "testing successful!" << endl;
}
