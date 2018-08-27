#include "nr3plus.h"

// if isign = 1, replaces data[0..2*n-1] by its discrete Fourier transform.
// if isign = -1, replaces data[0..2*n-1] by n times its inverse discrete Fourier transform.
// data is a complex array of length n stored as a real array of length 2*n.
// n must be an integer power of 2.
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