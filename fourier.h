#include "nr3plus.h"

void four1(Doub *data, Int_I n, Int_I isign);

void fft(VecComp_IO &data); 
void ifft(VecComp_IO &data);

void fft(MatComp_IO &data);
void ifft(MatComp_IO &data);

struct WrapVecDoub {
	VecDoub vvec;
	VecDoub &v;
	Int n, mask;

	WrapVecDoub(Int_I nn) : vvec(nn), v(vvec), n(nn/2),
	mask(n-1) {validate();}

	WrapVecDoub(VecDoub &vec) : v(vec), n(vec.size()/2),
	mask(n-1) {validate();}
		
	void validate() {if (n&(n-1)) throw("vec size must be power of 2");}

	inline Comp& operator[] (Int i) {return (Comp &)v[(i&mask) << 1];}

	inline Doub& real(Int i) {return v[(i&mask) << 1];}

	inline Doub& imag(Int i) {return v[((i&mask) << 1)+1];}

	operator VecDoub&() {return v;}
};

void realft(VecDoub_IO &data, Int_I isign);

void sinft(VecDoub_IO &y);

void cosft1(VecDoub_IO &y);

void cosft2(VecDoub_IO &y, Int_I isign);

Comp fft_interp(Doub_I x1, VecDoub_I &x, VecComp_I &y);

void fft_interp(VecComp_O &y1, VecDoub_I &x1, VecDoub_I &x, VecComp_I &y);

template <class T>
void fftshift(NRvector<T> &v)
{
	Long n{ v.size() };
	if (isodd(n)) error("fftshift only supports even columns!")
		Long i, halfn{ n / 2 };
	NRvector<T> temp(halfn);
	size_t size{ halfn * sizeof(T) };
	memcpy(&temp[0], v.ptr(), size);
	memcpy(v.ptr(), v.ptr() + halfn, size);
	memcpy(v.ptr() + halfn, &temp[0], size);
}

template <class T>
void fftshift(NRmatrix<T> &a, Int_I dim = 1)
{
	Long m{ a.nrows() }, n{ a.ncols() };
	if (dim == 1) {
		if (isodd(m)) error("fftshift only supports even rows!")
			Long halfm = m / 2;
		NRmatrix<T> temp(halfm, n);
		Long size = halfm*n * sizeof(T);
		memcpy(temp[0], a[0], size);
		memcpy(a[0], a[halfm], size);
		memcpy(a[halfm], temp[0], size);
	}
	else if (dim == 2) {
		if (isodd(n)) error("fftshift only supports even columns!")
			Long i, halfn{ n / 2 };
		NRvector<T> temp(halfn);
		Long size{ halfn * sizeof(T) };
		for (i = 0; i < m; ++i) {
			memcpy(&temp[0], a[i], size);
			memcpy(a[i], &a[i][halfn], size);
			memcpy(&a[i][halfn], &temp[0], size);
		}
	}
}

// === dft is strongly not recommanded, just read the sampling theorem and see how perfect fft is

// discrete fourier transform from X(x) to Y(k), no fftshift is needed
// each column of X is transformed to each column of Y
// using sum instead of integration, result not normalized
// for each column, Y_j = sum_i ( X_i*exp(-I*k_j*X_i) )
// this is much slower than fft, but for small (xmax-xmin) and (kmax-kmin), could be faster
void dft(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax);
void dft_par(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax);

// the inverse of dft, multiplied by 2*pi/(dx*dk).
// this might not be a precise inverse
void idft(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax);
void idft_par(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax);