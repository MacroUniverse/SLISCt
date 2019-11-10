// a wraper for cuBLAS
#include <cuda_runtime.h>
#include "cublas_v2.h"
#include <complex>
#include "../global.h"

namespace slisc {

// do zgemv with gpu (including allocation and data transfer)
void cuda_zgemv(Long M, Long N, const Comp *a, Long lda, const Comp *x, Long incx, Comp *y, Long incy)
{
	cublasHandle_t handle;
    cuDoubleComplex *devPtrA = 0, *devPtrX = 0, *devPtrY = 0;

    cudaMalloc((void**)&devPtrA, M*N*sizeof(Comp));
    cudaMalloc((void**)&devPtrX, N*sizeof(Comp));
    cudaMalloc((void**)&devPtrY, M*sizeof(Comp));

    cublasCreate(&handle);
    cublasSetMatrix(M, N, sizeof(Comp), a, M, devPtrA, M);
    cublasSetMatrix(N, 1, sizeof(Comp), x, M, devPtrX, M);
    cublasSetMatrix(M, 1, sizeof(Comp), y, M, devPtrY, M);

    Comp alpha(1, 0), beta(0, 0);
    cublasZgemv(handle, CUBLAS_OP_N, M, N, (cuDoubleComplex*)&alpha, devPtrA, M, devPtrX, incx, (cuDoubleComplex*)&beta, devPtrY, incy);

	cublasGetMatrix(M, 1, sizeof(Comp), devPtrY, M, y, M);
    cudaFree(devPtrA); cudaFree(devPtrX); cudaFree(devPtrY);
    cublasDestroy(handle);
}

} // namespace slisc
