#include "interp_2d.h"
#include "nr3plus.h"

void Spline_grid(MatComp_O &Z, VecDoub_I &x, VecDoub_I &y, MatComp_I &Z0, VecDoub_I &x0, VecDoub_I &y0)
{
	Long i, j, indi, indj, Nx0{ x0.size() }, Ny0{ y0.size() }, Nx{ x.size() }, Ny{ y.size() };
	Doub temp;
	if (Z.nrows() != Ny || Z.ncols() != Nx) Z.resize(Ny, Nx);
	// decide interpolation range
	const Doub x0min{ x0[0] }, x0max{ x0[Nx0 - 1] }, y0min{ y0[0] }, y0max{ y0[Ny0 - 1] };
	Long imin{ -1 }, Ni{ -1 }, jmin{ -1 }, Nj{ -1 };
	for (i = 0; i < Ny; ++i)
		if (y[i] >= y0min) {
			imin = i; break;
		}
	for (i = Ny - 1; i >= 0; --i)
		if (y[i] <= y0max) {
			Ni = i - imin + 1; break;
		}
	for (j = 0; j < Nx; ++j)
		if (x[j] >= x0min) {
			jmin = j; break;
		}
	for (j = Nx - 1; j >= 0; --j)
		if (x[j] <= x0max) {
			Nj = j - jmin + 1; break;
		}
	if (imin < 0 || Ni < 0 || jmin < 0 || Nj < 0) {
		Z = 0.; return;
	}
	MatComp Z1(Ny0, Nj); // for interp in x direction
	// interp. x direction
	{VecDoub Zx_R(Nx0), Zx_I(Nx0);
	for (i = 0; i < Ny0; ++i) {
		for (j = 0; j < Nx0; ++j) {
			Zx_R[j] = real(Z0[i][j]);
			Zx_I[j] = imag(Z0[i][j]);
		}
		Spline_interp f_R(x0, Zx_R), f_I(x0, Zx_I);
		for (j = 0; j < Nj; ++j) {
			temp = x[jmin + j];
			Z1[i][j] = Comp(f_R.interp(temp), f_I.interp(temp));
		}
	}}
	// interp. y direction
	{VecDoub Zy_R(Ny0), Zy_I(Ny0), y0flip, yflip;
	flip(y0flip, y0); flip(yflip, y);
	for (j = 0; j < Nj; ++j) {
		for (i = 0; i < Ny0; ++i) {
			Zy_R[i] = real(Z1[i][j]);
			Zy_I[i] = imag(Z1[i][j]);
		}
		Spline_interp f_R(y0flip, Zy_R), f_I(y0flip, Zy_I);
		indj = jmin + j;
		for (i = 0; i < Ni; ++i) {
			indi = imin + i;
			temp = yflip[indi];
			Z[indi][indj] = Comp(f_R.interp(temp), f_I.interp(temp));
		}
	}}
}
