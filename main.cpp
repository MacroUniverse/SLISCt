// comprehensive test of SLISC

//#include "test/test_all.h"

#include "SLISC/eig.h"
#include "SLISC/arithmetic.h"
#include "SLISC/disp.h"
//using namespace slisc;
//using namespace std;

namespace slisc {

// square diagonal matrix
// mostly a clone a Vector<T>
template <class T>
class Diag : public Vector<T>
{
private:
	typedef Vector<T> Base;
public:
	using Base::operator();
	using Base::operator=;
	Diag() : Base() {}
	Diag(Long_I N, const T &s) : Base(N, s) {}
	Diag(const Vector<T> &v) { *this = v; }
	Diag &operator=(const Diag &rhs)
	{ Base::operator=(rhs); return *this; }
	Diag &operator=(const Vector<T> &rhs)
	{ Base::operator=(rhs); return *this; }
	T &operator()(Long_I i, Long_I j)
	{
		if (i == j) return (*this)[i];
		return T();
	}
	const T &operator()(Long_I i, Long_I j) const
	{
		if (i == j) return (*this)[i];
		return T();
	}
};

// convert vector to diagonal matrix
template <class T>
const Diag<T> &diag(const Vector<T> &v)
{ return (Diag<T>&)v; }

template <class T, class T1, class T2>
void vec_scale(T *target, T1 *source, const T2 &s, Long_I N)
{
	for (Long i = 0; i < N; ++i) {
		target[i] = source[i] * s;
	}
}

// Cmat times Diag
template <class T, class T1, class T2>
void mul(Cmat<T> &v, const Cmat<T1> &v1, const Diag<T2> &v2)
{
	Long Nr = v1.nrows(), Nc = v1.ncols();
#ifdef _CHECKBOUNDS_
	if (Nc != v2.size()) error("illegal shape!");
#endif
	v.resize(Nr, v2.size());
	T * p = v.ptr();
	const T1 *p1 = v1.ptr();
	for (Long i = 0; i < Nc; ++i) {
		vec_scale(p, p1, v2[i], Nr);
		p += Nr; p1 += Nr;
	}
}

}

int main()
{
	using namespace slisc;
	// test_all();
	CmatDoub a(2, 2);
	a(0, 0) = 1.; a(1, 1) = 2.;
	a(1, 0) = 3.; a(0, 1) = 3.;
	CmatDoub eigVec; VecDoub eigVal;
	eigSym(eigVal, eigVec, a);
	disp(a, 10);
	disp(eigVec, 10);
	disp(eigVal, 10);
	CmatDoub eigVec1, eigVec2;
	mul(eigVec1, a, eigVec);
	disp(eigVec1, 10);
	mul(eigVec2, eigVec, diag(eigVal));
	disp(eigVec2, 10);
	eigVec1 -= eigVec2;
	std::cout << "max error = " << max(eigVec1) << std::endl;
}
