// test OpenMP
#pragma once
#include <omp.h>
#include "../SLISC/cmat.h"

// this program shows how to pre-allocate a workspace for each omp thread
void test_omp()
{
    using namespace slisc;
    printf("num_threas = %d\n", omp_get_num_threads());
    Int Nth = 3;
    Int N = 8;
    CmatInt wsp(N, Nth);
    CmatInt x1(N, N), x2(N, N);
    VecInt ys(N, 0), yp(N, 0); // serial and parallel solution
    for (Long i = 0; i < N * N; ++i) {
        x1[i] = randInt(30); x2[i] = randInt(30);
    }

    // serial solution
    // ys[j] = sum(x1(:,j) + x2(:,j));
    for (Long j = 0; j < N; ++j) {
        for (Long i = 0; i < N; ++i) {
            ys[j] += x1(i, j) + x2(i, j);
        }
    }
    
    // parallel solution using workspace for each thread
    omp_set_num_threads(Nth);
#pragma omp parallel for
    for (Long j = 0; j < N; ++j) {
        Int tid = omp_get_thread_num();
        printf("j = %lld thread = %d/%d\n", j, tid, omp_get_num_threads());
        for (Long i = 0; i < N; ++i) {
            wsp(i, tid) = x1(i, j) + x2(i, j);
        }
        yp[j] = sum(slice1(wsp, tid));
    }

    if (ys != yp)
        SLS_ERR("failed!");
    cout << "successful!" << endl;
}
