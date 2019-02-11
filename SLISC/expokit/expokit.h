#pragma once
#include "zgexpv.h"
#include "zhexpv.h"

namespace slisc {

// matrix / vector multiplication

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCoo<T> &a, Comp *x)
{
	mul_v_coo_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.size());
}

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCooH<T> &a, Comp *x)
{
	mul_v_cooh_v(y.ptr(), x.ptr(), a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.size());
}

// expv()

template <Char Option = 0, class T, SLS_IF(
	(is_MatCoo<T>() || is_MatCooH<T>()) &&
	(is_Comp<contain_type<T>>() || is_Doub<contain_type<T>>())
)>
inline void expv(VecComp_O v_out, const T &mat, VecComp_I v_in, Doub_I t, Int_I Nbase)
{
	const Doub tol = 0.; // set tolerance here
#ifdef SLS_CHECK_SHAPE
	if (mat.nrows() != mat.ncols() || mat.ncols() != v_in.size())
		error("wrong shape!");
#endif
	v_out.resize(v_in);
	Int iflag;
	VecComp wsp(MAX(10, SQR(mat.nrows()*(Nbase + 2) + 5 * (Nbase + 2)) + 7));
	VecInt iwsp(MAX(Nbase + 1, 7));

	if constexpr (Option == 'G' || Option == 0 && is_MatCoo<T>()) {
		ZGEXPV(v_in.size(), Nbase, t, v_in.ptr(), v_out.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
	else if (Option == 'H' || Option == 0 && is_MatCooH<T>()) {
		ZHEXPV(v_in.size(), Nbase, t, v_in.ptr(), v_out.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
}

} // namespace slisc
