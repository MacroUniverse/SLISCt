// classes for cuda matrix
#include "nr3.h"

Int Nbl, Nth; // <<<Nbl,Nth>>> for element-wise kernels

// always call this function first
void cuInit()
{
	Nbl = 320; Nth = 32;
}

// set elements to same value
// TODO use template <class T>
__global__
void cumemset(Doub* p, Long_I N, Doub s)
{
	Int i, ind, stride;
	ind = blockIdx.x * blockDim.x + threadIdx.x;
	stride = gridDim.x * blockDim.x;
	for (i = ind; i < N; i += stride)
		p[i] = s;
}

// device memory
template <class T>
class CUbase
{
protected:
	Long N;// number of elements
	T* p; // pointer to the first element
	size_t sizeT;
public:
	CUbase();
	CUbase(Long_I N);
	inline void resize(Long_I n);
	inline NRmat3d & operator=(const T &rhs);
	~CUbase();
}

template <class T>
CUbase<T>::CUbase(): N(0), p(nullptr), sizeT(sizeof(T)) {}

template <class T>
CUbase<T>::CUbase(Long_I n) : N(n), sizeT(sizeof(T))
{ cudaMalloc(&p, N*sizeT); }

template <class T>
inline void CUbase<T>::resize(Long_I n)
{
	if (n != N) {
		if (p != nullptr) cudaFree(p);
		N = n;
		if (n > 0)
			cudaMalloc(&p, N*sizeT);
		else
			p = nullptr;
	}
}

template <class T>
inline CUbase<T> & CUbase<T>::operator=(const T &rhs)
{
	if (N) cumemset<<<Nbl, Nth>>>(p, N, rhs);
	return *this;
}

template <class T>
NRbase<T>::~NRbase()
{ if (p) cudaFree(p); }
