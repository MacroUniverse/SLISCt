#pragma once
#include "expokit/zgexpv.h"
#include "expokit/zhexpv.h"

namespace slisc {

// matrix / vector multiplication

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCoo<T> &a, Comp *x)
{
	mul_v_coo_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.nnz());
}

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCooH<T> &a, Comp *x)
{
	mul_v_cooh_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.nrows(), a.nnz());
}

// expv()
// use ZGEXPV() for MatCoo<>, ZHEXPV() for MatCooH<>

template <Char Option = 0, class T, SLS_IF(
	(is_MatCoo<T>() || is_MatCooH<T>()) &&
	(is_Comp<contain_type<T>>() || is_Doub<contain_type<T>>())
)>
inline void expv(VecComp_IO v, const T &mat, Doub_I t, Int_I Nkrylov)
{
	const Doub tol = 0.; // set tolerance here
#ifdef SLS_CHECK_SHAPE
	if (mat.nrows() != mat.ncols() || mat.ncols() != v.size())
		error("wrong shape!");
#endif
	Int iflag;
	VecComp wsp(MAX(Long(10), SQR(mat.nrows()*(Nkrylov + 2) + 5 * (Nkrylov + 2)) + 7));
	VecInt iwsp(MAX(Nkrylov + 2, 7));

	if constexpr (Option == 'G' || Option == 0 && is_MatCoo<T>()) {
		ZGEXPV(v.size(), Nkrylov, t, v.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
	else if constexpr (Option == 'H' || Option == 0 && is_MatCooH<T>()) {
		ZHEXPV(v.size(), Nkrylov, t, v.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
	else error("unknown!");
}

} // namespace slisc
