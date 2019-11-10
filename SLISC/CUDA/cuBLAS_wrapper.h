#pragma once

namespace slisc {
	void cuda_zgemv(Long M, Long N, const Comp *a, Long lda, const Comp *x, Long incx, Comp *y, Long incy);
}
