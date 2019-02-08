#pragma once
#include "zgexpv.h"
#include "zhexpv.h"

namespace slisc {

template <class T>
inline void expGv(VecComp_O v_out, const MatCoo<T> &G, VecComp_I v_in, Doub_I t, Int_I Nbase)
{
	const tol = 0.; // set tolerance here
#ifdef SLS_CHECK_BOUNDS
	if (G.nrows() != G.ncols() || G.ncols() != v_in.size())
		error("wrong shape!");
#endif
	v_out.resize(v_in);
	Int iflag;
	VecComp wsp(MAX(10, SQR(n*(m + 2) + 5 * (m + 2)) + 7));
	VecInt iwsp(MAX(m + 1, 7));
	ZGEXPV(v_in.size(), Nbase, t, v_in.ptr(), v_out.ptr(),
		tol, norm_inf(G), wsp.ptr(), wsp.size(),
		iwsp.ptr(), iwsp.size(), G, 0, iflag);
}

template <class T>
inline void expHv(VecComp_O v_out, const MatCooH<T> &H, VecComp_I v_in, Doub_I t, Int_I Nbase)
{
	const tol = 0.; // set tolerance here
#ifdef SLS_CHECK_BOUNDS
	if (H.nrows() != H.ncols() || H.ncols() != v_in.size())
		error("wrong shape!");
#endif
	v_out.resize(v_in);
	Int iflag;
	VecComp wsp(MAX(10, SQR(n*(m + 2) + 5 * (m + 2)) + 7));
	VecInt iwsp(MAX(m + 1, 7));
	ZHEXPV(v_in.size(), Nbase, t, v_in.ptr(), v_out.ptr(),
		tol, norm_inf(H), wsp.ptr(), wsp.size(),
		iwsp.ptr(), iwsp.size(), H, 0, iflag);
}

}
