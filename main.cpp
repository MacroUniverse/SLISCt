// comprehensive test of nr3plus.h

#include "nr3.h"
#include "nr3plus.h"
#include "interp_1d.h"
#include "interp_2d.h"
#include "fourier.h"

using std::cout; using std::endl; using std::conj;


// test NRvector<>, NRmatrix<>, NRmat3d<> class themselves
void class_test()
{
	// default initialize
	{
	VecDoub vDoub;
	if (vDoub.size() != 0) error("failed!")
	if (vDoub.ptr() != nullptr) error("failed!")
	MatDoub aDoub;
	if (aDoub.size() != 0) error("failed!")
	if (aDoub.nrows() != 0) error("failed!")
	if (aDoub.ncols() != 0) error("failed!")
	if (aDoub.ptr() != nullptr) error("failed!")
	Mat3Doub a3Doub;
	if (a3Doub.size() != 0) error("failed!")
	if (a3Doub.dim1() != 0) error("failed!")
	if (a3Doub.dim2() != 0) error("failed!")
	if (a3Doub.dim3() != 0) error("failed!")
	if (a3Doub.ptr() != nullptr) error("failed!")
	}

	// size initialize
	{
	VecDoub vDoub(3);
	if (vDoub.size() != 3) error("failed!")
	if (vDoub.ptr() != &vDoub[0]) error("failed!")
	MatDoub aDoub(3, 3);
	if (aDoub.size() != 9) error("failed!")
	if (aDoub.nrows() != 3) error("failed!")
	if (aDoub.ncols() != 3) error("failed!")
	if (aDoub.ptr() != &aDoub[0][0]) error("failed!")
	Mat3Doub a3Doub(3, 3, 3);
	if (a3Doub.size() != 27) error("failed!")
	if (a3Doub.dim1() != 3) error("failed!")
	if (a3Doub.dim2() != 3) error("failed!")
	if (a3Doub.dim3() != 3) error("failed!")
	if (a3Doub.ptr() != &a3Doub[0][0][0]) error("failed!")
	}

	// const initialize
	VecDoub vDoub(3, 1.);
	if (vDoub != 1.) error("failed!")
	MatDoub aDoub(3, 3, 1.);
	if (aDoub != 1.) error("failed!")
	Mat3Doub a3Doub(3, 3, 3, 1.);
	if (a3Doub != 1.) error("failed!")

	// resize
	vDoub.resize(0);
	if (vDoub.ptr() != nullptr) error("failed!")
	vDoub.resize(4);
	if (vDoub.size() != 4) error("failed!")
	if (vDoub.ptr() != &vDoub[0]) error("failed!")
	aDoub.resize(0, 3);
	if (aDoub.ptr() != nullptr) error("failed!")
	aDoub.resize(3, 0);
	if (aDoub.ptr() != nullptr) error("failed!")
	aDoub.resize(4, 4);
	if (aDoub.size() != 16) error("failed!")
	if (aDoub.nrows() != 4) error("failed!")
	if (aDoub.ncols() != 4) error("failed!")
	if (aDoub.ptr() != &aDoub[0][0]) error("failed!")
	a3Doub.resize(0, 0, 4);
	if (a3Doub.ptr() != nullptr) error("failed!")
	a3Doub.resize(0, 4, 0);
	if (a3Doub.ptr() != nullptr) error("failed!")
	a3Doub.resize(4, 0, 0);
	if (a3Doub.ptr() != nullptr) error("failed!")
	a3Doub.resize(4, 4, 4);
	if (a3Doub.size() != 64) error("failed!")
	if (a3Doub.dim1() != 4) error("failed!")
	if (a3Doub.dim2() != 4) error("failed!")
	if (a3Doub.dim3() != 4) error("failed!")
	if (a3Doub.ptr() != &a3Doub[0][0][0]) error("failed!")

	// assignment operator
	vDoub = 1.; if (vDoub != 1.) error("failed!")
	aDoub = 1.; if (aDoub != 1.) error("failed!")
	a3Doub = 1.; if (a3Doub != 1.) error("failed!")
	VecDoub vDoub1(4);
	vDoub1 = 2.;
	vDoub = vDoub1;
	if (vDoub != vDoub1) error("failed!")
	MatDoub aDoub1(4, 4);
	aDoub1 = 2.;
	aDoub = aDoub1;
	if (aDoub != aDoub1) error("failed!")
	Mat3Doub a3Doub1(4, 4, 4);
	a3Doub1 = 2.;
	a3Doub = a3Doub1;
	if (a3Doub != a3Doub1) error("failed!")

	// move operator
	VecDoub vDoub2;
	MatDoub aDoub2;
	Mat3Doub a3Doub2;
	vDoub2 << vDoub;
	if (vDoub2 != vDoub1) error("failed!")
	if (vDoub.size() != 0) error("failed!")
	if (vDoub.ptr() != 0)  error("failed!")
	aDoub2 << aDoub;
	if (aDoub2 != aDoub1)  error("failed!")
	if (aDoub.size() != 0) error("failed!")
	if (aDoub.nrows() != 0) error("failed!")
	if (aDoub.ncols() != 0) error("failed!")
	a3Doub2 << a3Doub;
	if (a3Doub2 != a3Doub1) error("failed!")
	if (a3Doub.size() != 0) error("failed!")
	if (a3Doub.dim1() != 0) error("failed!")
	if (a3Doub.dim2() != 0) error("failed!")
	if (a3Doub.dim3() != 0) error("failed!")

	// end()
	vDoub1.end() = 3.1;
	if (vDoub1.end() != 3.1)  error("failed!")
	aDoub1.end() = 3.1;
	if (aDoub1.end() != 3.1)  error("failed!")
	a3Doub1.end() = 3.1;
	if (a3Doub1.end() != 3.1)  error("failed!")

	// operator[]
	vDoub1[vDoub1.size()-1] = 5.5;
	if ( vDoub1[vDoub1.size()-1] != 5.5 ) error("failed!")
	if (vDoub1.end() != 5.5)  error("failed!")
	aDoub1[aDoub1.nrows()-1][aDoub1.ncols()-1] = 5.5;
	if ( aDoub1[aDoub1.nrows()-1][aDoub1.ncols()-1] != 5.5 ) error("failed!")
	if (aDoub1.end() != 5.5)  error("failed!")
	a3Doub1[a3Doub1.dim1()-1][a3Doub1.dim2()-1][a3Doub1.dim3()-1] = 5.5;
	if ( a3Doub1[a3Doub1.dim1()-1][a3Doub1.dim2()-1][a3Doub1.dim3()-1] != 5.5 ) error("failed!")
	if (a3Doub1.end() != 5.5)  error("failed!")

	// operator()
	// TODO:
}

// test +=, -=, *=, /=
void test_self_op()
{
	VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
	VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
	VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

	// vector ?= scalar

	vLlong = 1;
	vLlong += 1;
	if (vLlong != 2) error("failed!")
	vLlong -= 1;
	if (vLlong != 1) error("failed!")
	vLlong *= 2;
	if (vLlong != 2) error("failed!")
	// TODO: vLlong /= 2; failed!

	vDoub = 1.;
	vDoub += 1.;
	if (vDoub != 2.) error("failed!")
	vDoub -= 1.;
	if (vDoub != 1.) error("failed!")
	vDoub *= 2.;
	if (vDoub != 2.) error("failed!")
	vDoub /= 2.;
	if (vDoub != 1.) error("failed!")

	vComp = Comp(1., 1.);
	vComp += Comp(1., 1.);
	if (vComp != Comp(2., 2.)) error("failed!")
	vComp -= Comp(1., 1.);
	if (vComp != Comp(1., 1.)) error("failed!")
	vComp *= 2.;
	if (vComp != Comp(2., 2.)) error("failed!")
	vComp /= 2.;
	if (vComp != Comp(1., 1.)) error("failed!")

	// vector ?= vector
	
	vLlong = 1; vLlong1 = 1;
	vLlong += vLlong1;
	if (vLlong != 2) error("failed!")
	vLlong -= vLlong1;
	if (vLlong != 1) error("failed!")
	vLlong2 = 2;
	vLlong *= vLlong2;
	if (vLlong != 2) error("failed!")
	vLlong /= vLlong2;
	if (vLlong != 1) error("failed!")

	vDoub = 1.; vDoub1 = 1.;
	vDoub += vDoub1;
	if (vDoub != 2.) error("failed!")
	vDoub -= vDoub1;
	if (vDoub != 1.) error("failed!")
	vDoub2 = 2.;
	vDoub *= vDoub2;
	if (vDoub != 2.) error("failed!")
	vDoub /= vDoub2;
	if (vDoub != 1.) error("failed!")

	vComp = Comp(1., 1.); vComp1 = Comp(1., 1.);
	vComp += vComp1;
	if (vComp != Comp(2., 2.)) error("failed!")
	vComp -= vComp1;
	if (vComp != Comp(1., 1.)) error("failed!")
	vComp2 = 2.;
	vComp *= vComp2;
	if (vComp != Comp(2., 2.)) error("failed!")
	vComp /= vComp2;
	if (vComp != Comp(1., 1.)) error("failed!")
}

void test_plus_minus_times_devide()
{
	VecInt vLlong(3), vLlong1(3), vLlong2(3), vLlong3(3);
	VecDoub vDoub(3), vDoub1(3), vDoub2(3), vDoub3(3);
	VecComp vComp(3), vComp1(3), vComp2(3), vComp3(3);

	// op(v, v, s)

	// TODO: Llong version
	vDoub1 = 1.;
	plus(vDoub, vDoub1, 1.);
	if (vDoub != 2.) error("failed!")
	minus(vDoub, vDoub1, 1.);
	if (vDoub != 0.) error("failed!")
	times(vDoub, vDoub1, 2.);
	if (vDoub != 2.) error("failed!")
	divide(vDoub, vDoub1, 2.);
	if (vDoub != 0.5) error("failed!")
	// TODO: Comp version
	vDoub = 1.;
	plus(vComp, vDoub1, Comp(1., 1.));
	if (vComp != Comp(2., 1.)) error("failed!")

	// op(v, s, v)

	// TODO: Llong version
	vDoub1 = 1.;
	plus(vDoub, 1., vDoub1);
	if (vDoub != 2.) error("failed!")
	minus(vDoub, 1., vDoub1);
	if (vDoub != 0.) error("failed!")
	times(vDoub, 2., vDoub1);
	if (vDoub != 2.) error("failed!")
	vDoub1 = 2.;
	divide(vDoub, 2., vDoub1);
	if (vDoub != 1.) error("failed!")
	// TODO: Comp version
	vDoub1 = 1.;
	plus(vComp, vDoub1, Comp(1., 1.));
	if (vComp != Comp(2., 1.)) error("failed!")

	// op(v, v, v)

	vLlong1 = 4; vLlong2 = 2;
	plus(vLlong, vLlong1, vLlong2);
	if (vLlong != 6) error("failed!")
	minus(vLlong, vLlong1, vLlong2);
	if (vLlong != 2) error("failed!")
	times(vLlong, vLlong1, vLlong2);
	if (vLlong != 8) error("failed!")
	divide(vLlong, vLlong1, vLlong2);
	if (vLlong != 2) error("failed!")

	vDoub1 = 1.; vDoub2 = 2.;
	plus(vDoub, vDoub1, vDoub2);
	if (vDoub != 3.) error("failed!")
	minus(vDoub, vDoub2, vDoub1);
	if (vDoub != 1.) error("failed!")
	times(vDoub, vDoub1, vDoub2);
	if (vDoub != 2.) error("failed!")
	divide(vDoub, vDoub1, vDoub2);
	if (vDoub != 0.5) error("failed!")
	
	vComp1 = Comp(1., 1.); vComp2 = Comp(2., 2.);
	plus(vComp, vComp1, vComp2);
	if (vComp != Comp(3., 3.)) error("failed!")
	minus(vComp, vComp2, vComp1);
	if (vComp != Comp(1., 1.)) error("failed!")
	vComp2 = 2.;
	times(vComp, vComp1, vComp2);
	if (vComp != Comp(2., 2.)) error("failed!")
	divide(vComp, vComp1, vComp2);
	if (vComp != Comp(0.5, 0.5)) error("failed!")
}

void fft_test()
{
	// fft(VecComp_IO) and ifft(VecComp_IO)
	VecComp v(4); v[0] = 1; v[1] = I;  v[2] = -1; v[3] = -I;
	fft(v);
	VecComp v1(4, 0.); v1[1] = 4.;
	v -= v1;
	if (max(v) > 1.5e-15) error("failed!")
	ifft(v1);
	v[0] = 1; v[1] = I;  v[2] = -1; v[3] = -I;
	v1 /= 4.; v1 -= v;
	if (max(v1) > 1e-15) error("failed!")


	// fft_interp()
	VecDoub x; linspace(x, 1., 3., 3);
	VecComp y; linspace(y, Comp(1., 1.), Comp(3., 3.), 3);
	if (abs(fft_interp(x[0], x, y) - y[0]) > 1e-15) error("failed!");
	if (abs(fft_interp(x[1], x, y) - y[1]) > 1e-15) error("failed!");
	if (abs(fft_interp(x[2], x, y) - y[2]) > 1e-15) error("failed!");

	// fftshift()
	VecInt vInt; linspace(vInt, 1, 4, 4);
	fftshift(vInt);
	VecInt vInt1(4); vInt1[0] = 3; vInt1[1] = 4; vInt1[2] = 1; vInt1[3] = 2;
	if (vInt != vInt1) error("failed!");
}

void four2x(Doub *data2, Doub *data, Int_I n, Int_I isign) {
	Int nn,mmax,m,j,istep,i;
	Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) error("n must be power of 2 in four1")
	nn = n << 1;
	// get bit inverse order
	j = 1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			SWAP(data[j-1],data[i-1]);
			SWAP(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	// do two element dft
	for (i = 0; i < 2*n; i += 4) {
		j = 2*i;
		data2[j+2] = data2[j] = data[i];
		data2[j+3] = data2[j+1] = data[i+1];
		data2[j+6] = -(data2[j+4] = data[i+2]);
		data2[j+7] = -(data2[j+5] = data[i+3]);
	}

	// do the rest
	nn <<= 1;
	mmax=4;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data2[j-1]-wi*data2[j];
				tempi=wr*data2[j]+wi*data2[j-1];
				data2[j-1]=data2[i-1]-tempr;
				data2[j]=data2[i]-tempi;
				data2[i-1] += tempr;
				data2[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fft2x(VecComp_O &data2, VecComp_I &data)
{
	VecComp temp; temp = data;
	data2.resize(data.size()*2);
	four2x((Doub*)data2.ptr(), (Doub*)temp.ptr(), data.size(), -1);
}

void four4x(Doub *data2, Doub *data, Int_I n, Int_I isign) {
	Int nn,mmax,m,j,istep,i;
	Doub wtemp,wr,wpr,wpi,wi,theta,tempr,tempi;
	if (n<2 || n&(n-1)) error("n must be power of 2 in four1")
	nn = n << 1;
	// get bit inverse order
	j = 1;
	for (i=1;i<nn;i+=2) {
		if (j > i) {
			SWAP(data[j-1],data[i-1]);
			SWAP(data[j],data[i]);
		}
		m=n;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	// do two element dft
	for (i = 0; i < 4*n; i += 8) {
		j = 4*i;
		data2[j+6] = data2[j+4] = data2[j+2] = data2[j] = data[i];
		data2[j+7] = data2[j+5] = data2[j+3] = data2[j+1] = data[i+1];
		data2[j+6] = -(data2[j+4] = data[i+2]);
		data2[j+7] = -(data2[j+5] = data[i+3]);
		data2[j+12] = data2[j+11] = -(data2[j+15] = data2[j+8] = data[i+2]);
		data2[j+14] = data2[j+13] = -(data2[j+10] = data2[j+9] = data[i+3]);
	}

	// do the rest
	nn <<= 2;
	mmax=8;
	while (nn > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=nn;i+=istep) {
				j=i+mmax;
				tempr=wr*data2[j-1]-wi*data2[j];
				tempi=wr*data2[j]+wi*data2[j-1];
				data2[j-1]=data2[i-1]-tempr;
				data2[j]=data2[i]-tempi;
				data2[i-1] += tempr;
				data2[i] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}

void fft4x(VecComp_O &data2, VecComp_I &data)
{
	VecComp temp; temp = data;
	data2.resize(data.size()*4);
	four2x((Doub*)data2.ptr(), (Doub*)temp.ptr(), data.size(), -1);
}

// new test scratch
void test()
{
	VecComp v(4); v[0] = Comp(1., 1.); v[1] = Comp(2., 2.); v[2] = Comp(3., 5.); v[3] = Comp(4., 7.);
	disp(v);
	VecComp v1; v1 = v;
	fft(v1);
	disp(v1);
	v1 = v;
	VecComp v2;
	fft2x(v2, v1);
	disp(v2);
	VecComp v3(8, 0.); v3[0] = v[0]; v3[1] = v[1]; v3[6] = v[2]; v3[7] = v[3];
	fft(v3);
	disp(v3);
}

int main()
{
	// temporary test
	test();

	// systematic tests
	test_self_op();
	test_plus_minus_times_devide();
	class_test();
	fft_test();

	// test operator() and end()
	VecDoub xv, xv1(3), xv2(3);
	linspace(xv1, 1., 3.); linspace(xv2, 4., 6.);

	//// test new disp()
	//VecUchar v8(3);
	//linspace<Uchar>(v8, 1, 3);
	//VecInt vi(3);
	//linspace<Int>(vi, 1, 3);
	//MatUchar A8(2, 3);
	//linspace<Uchar>(A8, 1, 6);
	//MatInt AI(2, 3);
	//linspace<Int>(AI, 1, 6);
	//Mat3Doub A3(2, 2, 2);
	//linspace<Doub>(A3, 1., 8.);
	//Mat3Comp C3;
	//C3.resize(2, 2, 2);
	//linspace<Comp>(C3, Comp(0, 1), Comp(14, 15));
	//
	//// test cubic spline interp 1D
	//Int N{ 2 };
	//VecDoub xx(N), yy(N);
	//linspace(xx, 0., 1.); linspace(yy, 1., 2.);
	//Spline_interp myfun(xx, yy);
	//cout << myfun.interp(0.5) << endl;

	//// test cubic spline interp 2D
	//Int i, j, m{ 5 }, n{ 5 };
	//MatComp Psi(m, n), Psi0(m, n);
	//VecDoub x0, y0, x, y;
	//linspace<Doub>(x0, 1, n, n); linspace<Doub>(y0, 1, m, m);
	//for (i = 0; i < m; ++i)
	//	for (j = 0; j < n; ++j)
	//		Psi0[i][j] = Comp(x0[j]*x0[j], y0[n - i - 1]*y0[n - i - 1]);
	//disp(Psi0);
	//x = x0; y = y0;
	//Spline_grid(Psi, x, y, Psi0, x0, y0);
	//disp(Psi);

	//VecDoub v(3,0.); VecComp vc(3,0.);
	//linspace(v, 0., PI);
	////disp(v, 8);
	////disp(v, 0, 2, 8);
	//linspace<Comp>(vc, 0., Comp(PI, 2.*PI));
	////disp(vc, 8);
	////disp(vc, 0, 2, 8);
	//Comp c, c1(3.,5.);
	//c = c1 + 1;
	//cout << c << endl;
	//c = c1 - 1;
	//cout << c << endl;
	//c = 1 - c1;
	//cout << c << endl;
	//c = c1 * 2;
	//cout << c << endl;
	//c = 2 * c1;
	//cout << c << endl;
	//c = 2/c1;
	//cout << c << endl;
	//MatDoub a(3,3,0.);
	//linspace(a, 0., 8.);
	////disp(a, 5);
}
