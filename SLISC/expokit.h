#pragma once
#include "expokit/zgexpv.h"
#include "expokit/zhexpv.h"

namespace slisc {

// matrix / vector multiplication

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCoo<T> &a, Comp *x)
{
	mul_v_coo_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

template <class T, SLS_IF(is_scalar<T>())>
void mul(Comp *y, const MatCooH<T> &a, Comp *x)
{
	mul_v_cooh_v(y, x, a.ptr(), a.row_ptr(), a.col_ptr(), a.n1(), a.nnz());
}

// expv()
// use ZGEXPV() for MatCoo<>, ZHEXPV() for MatCooH<>

template <Char Option = 0, class Tvec, class Tmat, SLS_IF(
	is_dense_vec<Tvec>() &&
	(is_Comp<contain_type<Tvec>>() || is_Doub<contain_type<Tvec>>()) &&
	(is_MatCoo<Tmat>() || is_MatCooH<Tmat>()) &&
	(is_Comp<contain_type<Tmat>>() || is_Doub<contain_type<Tmat>>())
)>
inline void expv(Tvec &v, const Tmat &mat, Doub_I t, Int_I Nkrylov)
{
	const Doub tol = 0.; // set tolerance here
#ifdef SLS_CHECK_SHAPE
	if (mat.n1() != mat.n2() || mat.n2() != v.size())
		SLS_ERR("wrong shape!");
#endif
	Int iflag;
	VecComp wsp(MAX(Long(10), SQR(mat.n1()*(Nkrylov + 2) + 5 * (Nkrylov + 2)) + 7));
	VecInt iwsp(MAX(Nkrylov + 2, 7));

	if constexpr (Option == 'G' || Option == 0 && is_MatCoo<Tmat>()) {
		ZGEXPV(v.size(), Nkrylov, t, v.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
	else if constexpr (Option == 'H' || Option == 0 && is_MatCooH<Tmat>()) {
		ZHEXPV(v.size(), Nkrylov, t, v.ptr(),
			tol, norm_inf(mat), wsp.ptr(), wsp.size(),
			iwsp.ptr(), iwsp.size(), mat, 0, iflag);
	}
	else SLS_ERR("unknown!");
}

} // namespace slisc
