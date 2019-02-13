#pragma once
#include "cmat.h"
#include "sparse.h"

namespace slisc {

// generate Gauss-Lobatto abscissas and weights
// pt is in [-1,1]
// data from https://keisan.casio.com/exec/system/1280801905
// use 38 digits data in case long double is needed
inline void GaussLobatto(VecDoub_O &x, VecDoub_O &w, Int_I N)
{
	Int i, N2 = N/2;
	x.resize(N); w.resize(N);
	x[0] = -1; x[N - 1] = 1;
	if (N == 10) {
		x[5] = 0.16527895766638702462621976595817353323;
		x[6] = 0.47792494981044449566117509273125799789;
		x[7] = 0.73877386510550507500310617485983072502;
		x[8] = 0.91953390816645881382893266082233813415;

		w[5] = 0.32753976118389745665651052791689314471;
		w[6] = 0.29204268367968375787558225737444389221;
		w[7] = 0.22488934206312645211945782173104784275;
		w[8] = 0.13330599085107011112622717075539289811;
		w[9] = 0.022222222222222222222222222222222222222;
	}
	else {
		printf("no data!"); exit(EXIT_FAILURE);
	}
	for (i = 1; i < N2; ++i)
		x[i] = -x[N - i - 1];
	for (i = 0; i < N2; ++i)
		w[i] = w[N - i - 1];
}


// get derivatives of legendre interpolation polynomials t abscissas
// df is an NxN matrix
// f'_i(x_j) = df(j,i)
inline void legendre_interp_der(CmatDoub_O &df, VecDoub_I &x)
{
	Long i, j, k, N{ x.size() };
	Doub t;
	df.resize(N,N);
	for (i = 0; i < N; ++i)
		for (j = 0; j < N; ++j) {
			if (j != i) {
				t = 1.;
				for (k = 0; k < i; ++k) {
					if (k == j) t /= x[i] - x[k];
					else t *= (x[j] - x[k]) / (x[i] - x[k]);
				}
				for (k = i + 1; k < N; ++k) {
					if (k == j) t /= x[i] - x[k];
					else t *= (x[j] - x[k]) / (x[i] - x[k]);
				}
			}
			else {
				t = 0.;
				for (k = 0; k < i; ++k)
					t += 1. / (x[i] - x[k]);
				for (k = i + 1; k < N; ++k)
					t += 1. / (x[i] - x[k]);
			}
			df(j,i) = t;
		}
}

// calculate FEDVR global index from FE index i and DVR index j
inline Long indFEDVR(Long_I i, Long_I j, Long_I Ngs)
{ return (Ngs-1) * i + j - 1; }

// generate FEDVR grid and weight
// 'Nfe' is the number of finite elements
// 'Ngs' is the number of abscissas in each element (including 2 boundaries)
// 'bounds' size = Nfe+1, are the boundary points of finite elements
// 'x0','w0' are grid points and weights in [-1,1]
// 'x','w' are the global grid points and weights

inline void FEDVR_grid(VecDoub_O &x, VecDoub_O &w, VecDoub_I &wFE, VecDoub_I &xFE, VecDoub_I &x0, VecDoub_I &w0)
{
	Long i, j, k = 0;
	Long Ngs = x0.size();
	Long Nfe = wFE.size();
	Long Nx = Nfe * (Ngs - 1) - 1;
	Doub a, b;
	x.resize(Nx); w.resize(Nx);
	for (i = 0; i < Nfe; ++i) {
		a = wFE(i);
		b = xFE(i);
		if (i > 0)
			w[k-1] += a*w0(0);
		for (j = 1; j < Ngs; ++j) {
			x[k] = a*x0[j] + b;
			w[k] = a*w0[j];
			if (++k == Nx) break;
		}
	}
}

// generate dense kinetic matrix
inline void D2_matrix(McooDoub_O &D2, VecDoub_I &w0, VecDoub_I &wFE, CmatDoub_I &df)
{
	Long i, j, k, m, n, mm, nn;
	Long Nfe = wFE.size();
	Long Ngs = w0.size();
	Long N = Nfe * (Ngs - 1) - 1;
	Doub s, coeff;

	// prepare basic block
	CmatDoub block(Ngs, Ngs);
	for (j = 0; j < Ngs; ++j) {
		for (i = 0; i <= j; ++i) {
			s = 0.;
			for (k = 0; k < Ngs; ++k)
				s += w0(k) * df(k, i) * df(k, j);
			block(i, j) = s;
		}
	}

	// calculate Kinetic matrix T
	D2.reshape(N, N);
	D2.reserve((Ngs*Ngs - 1)*Nfe - 4 * Ngs + 3); // # non-zero elements
	for (i = 0; i < Nfe; ++i) {
		// blocks without boundary
		for (n = 1; n < Ngs - 1; ++n) {
			coeff = -1. / SQR(wFE(i));
			for (m = 1; m <= n; ++m) {
				mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i, n, Ngs);
				s = coeff * block(m,n);
				D2.push(s, mm, nn);
				if (mm != nn) {
					D2.push(s, nn, mm);
				}
			}
		}
	}

	for (i = 0; i < Nfe - 1; ++i) {
		// block right boundary
		n = Ngs - 1;
		coeff = -1. / (pow(wFE(i), 1.5) * ::sqrt(wFE(i) + wFE(i + 1)));
		for (m = 1; m < n; ++m) {
			mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i, n, Ngs);
			s = coeff * block(m, n);
			D2.push(s, mm, nn); D2.push(s, nn, mm);
		}

		// block lower right corner
		m = Ngs - 1;
		coeff = -1. / (wFE(i) + wFE(i + 1));
		mm = indFEDVR(i, m, Ngs);
		s = coeff*(block(m,n)/wFE(i) + block(0,0)/wFE(i+1));
		D2.push(s, mm, mm);

		// block upper boundary
		coeff = -1. / (pow(wFE(i + 1), 1.5) * ::sqrt(wFE(i) + wFE(i + 1)));
		for (n = 1; n < Ngs - 1; ++n) {
			mm = indFEDVR(i, m, Ngs); nn = indFEDVR(i + 1, n, Ngs);
			s = coeff * block(0, n);
			D2.push(s, mm, nn); D2.push(s, nn, mm);
		}
	}

	for (i = 0; i < Nfe - 2; ++i) {
		// block upper right corner
		coeff = -1. / (wFE(i + 1) * ::sqrt((wFE(i) + wFE(i + 1))*(wFE(i + 1) + wFE(i + 2))));
		mm = indFEDVR(i, Ngs-1, Ngs); nn = indFEDVR(i + 1, Ngs-1, Ngs);
		s = coeff * block(0, n);
		D2.push(s, mm, nn); D2.push(s, nn, mm);
	}
}

// bounds: FE boundaries, size = Nfe + 1
// Ngs: grid points per finite element (including boundaries)
void D2_matrix(McooDoub_O D2, VecDoub_O x, VecDoub_O w, VecDoub_O u, VecDoub_I bounds, Int_I Ngs)
{
	Int i, j;
	Int Nfe = bounds.size() - 1; // number of finite elements
	Int Nx = (Ngs - 1)*Nfe - 1; // total grid points

	// grid points, weights, base function values in [-1, 1]
	VecDoub x0, w0, f0;
	GaussLobatto(x0, w0, Ngs);
	invSqrt(f0, w0);
	CmatDoub df(Ngs, Ngs); // df(i, j) = f_j(x_i)
	legendre_interp_der(df, x0);
	for (i = 0; i < Ngs; ++i)
		for (j = 0; j < Ngs; ++j)
			df(j, i) *= f0(i);

	// midpoints, widths and bounds of finite elements
	VecDoub xFE(Nfe), wFE(Nfe);

	for (i = 0; i < Nfe; ++i) {
		wFE(i) = 0.5*(bounds(i + 1) - bounds(i));
		xFE(i) = 0.5*(bounds(i) + bounds(i + 1));
	}

	FEDVR_grid(x, w, wFE, xFE, x0, w0);
	invSqrt(u, w);

	// Sparse Hamiltonian
	D2_matrix(D2, w0, wFE, df);
}
} // namespace slisc
