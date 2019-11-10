#include "../time.h"
#include "../arithmetic.h"
#include "../cmat.h"
#include "../disp.h"
#include "cuBLAS_wrapper.h"

using namespace slisc;

int main()
{
	Long M = 4, N = 3;
	CmatComp a(M, N);
	VecComp x(N), y(M), y1(M);
	linspace(x, Comp(1, 2), Comp(N, N+1));
	linspace(a, Comp(1, 2), Comp(M*N, M*N+1));
	cuda_zgemv(M, N, a.ptr(), M, x.ptr(), 1, y.ptr(), 1);

	mul(y1, a, x);
	y1 -= y;
	if (max_abs(y1) > 0)
		SLS_ERR("failed!");
	cout << "test successful!" << endl;
}
