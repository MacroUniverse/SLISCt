// standard read write
#include "benchmark.h"
using std::cout; using std::endl; using std::conj;

void bench_read_write()
{
	Long i, N1 = 10000, N2 = 10000, N = N1 * N2;
	MatDoub A(N1, N2), B(N1, N2);
	Doub *pa = A.ptr(), *pb = B.ptr();
	ctic();
	for (i = 0; i < N; ++i) {
		pa[i] = (Doub)N;
	}
	cout << ctoc() << endl;
}
