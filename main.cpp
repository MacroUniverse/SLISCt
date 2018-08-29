// comprehensive test of nr3plus.h

#include "nr3.h"
#include "nr3plus.h"
#include "interp_1d.h"
#include "interp_2d.h"
#include "fourier.h"

using std::cout; using std::endl; using std::conj;


// test NRvector<>, NRmatrix<>, NRmat3d<> class themselves
void test_class()
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
	vDoub1[vDoub1.size()-1] = 3.1;
	if (vDoub1.end() != 3.1)  error("failed!")
	if (vDoub1.end(1) != 3.1)  error("failed!")
	aDoub1(aDoub1.size()-1) = 3.1;
	if (aDoub1.end() != 3.1)  error("failed!")
	if (aDoub1.end(1) != 3.1)  error("failed!")
	a3Doub1(a3Doub1.size()-1) = 3.1;
	if (a3Doub1.end() != 3.1)  error("failed!")
	if (a3Doub1.end(1) != 3.1)  error("failed!")
	vDoub1[vDoub1.size()-2] = 3.2;
	if (vDoub1.end(2) != 3.2)  error("failed!")
	aDoub1(aDoub1.size()-2) = 3.2;
	if (aDoub1.end(2) != 3.2)  error("failed!")
	a3Doub1(a3Doub1.size()-2) = 3.2;
	if (a3Doub1.end(2) != 3.2)  error("failed!")

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

// test plus(), minus(), times(), devide()
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

// test fft module
void test_fft()
{
	// test bit_inv()
	VecComp v; linspace(v, 1., 16., 16);
	VecComp v1(16);
	bit_inv(v1.ptr(), v.ptr(), v.size());
	VecComp v2; v2 = v;
	bit_inv(v2.ptr(), v2.size());
	if (v1 != v2) error("failed!")
	bit_inv(v1.ptr(), v1.size());
	if (v1 != v) error("failed!")

	// fft(VecComp_IO) and ifft(VecComp_IO)
	v.resize(4); v[0] = 1; v[1] = I;  v[2] = -1; v[3] = -I;
	fft(v);
	v1.resize(4); v1 = 0.; v1[1] = 4.;
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

	// test fft2x(), ifft2x, fft4x(), ifft4x
	v.resize(4); v[0] = Comp(1., 1.); v[1] = Comp(2., 2.); v[2] = Comp(3., 5.); v[3] = Comp(4., 7.);
	fft2x(v2, v);
	VecComp v3(8, 0.); v3[0] = v[0]; v3[1] = v[1]; v3[6] = v[2]; v3[7] = v[3];
	fft(v3);
	if (v2 != v3) error("failed!")

	ifft2x(v2, v);
	v3 = 0.; v3[0] = v[0]; v3[1] = v[1]; v3[6] = v[2]; v3[7] = v[3];
	ifft(v3);
	if (v2 != v3) error("failed!")

	VecComp v4;
	fft4x(v4, v);
	VecComp v5(16, 0.);  v5[0] = v[0]; v5[1] = v[1]; v5[14] = v[2]; v5[15] = v[3];
	fft(v5);
	v4 -= v5;
	if (max(v4) > 1e-14) error("failed!")

	ifft4x(v4, v);
	v5 = 0.; v5[0] = v[0]; v5[1] = v[1]; v5[14] = v[2]; v5[15] = v[3];
	ifft(v5);
	v4 -= v5;
	if (max(v4) > 1e-14) error("failed!")
}

// new test scratch
void test()
{
	
}

int main()
{
	// temporary test
	test();

	// systematic tests
	test_class();
	test_self_op();
	test_plus_minus_times_devide();
	test_fft();

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
