// include all header files
#pragma once

// === all macros ===
// SLS_CPP_STD = 17
// SLS_CHECK_BOUNDS
// SLS_CHECK_SHAPE
// SLS_RAND_SEED
// SLS_USE_INT_AS_LONG (this one must be defined for every translation unit, i.e. in "global.h" or in compiler option)
// SLS_CHECK_COO_REPEAT
// SLS_ALLOW_COPY_CONSTRUCTOR
// SLS_CUSLISC
// SLS_USE_MKL
// SLS_USE_CBLAS
// SLS_USE_LAPACKE
// SLS_USE_GSL
// SLS_FP_EXCEPT
// SLS_USE_UTFCPP
// SLS_ERR
// SLS_WARN
// SLS_IF0
// SLS_IF
// SLS_TIME_H_ERR
// SLS_MATT_REPLACE
// SLS_HAS_FILESYSTEM (only define when <filesystem> works, don't define for linux)

// basics
#include "global.h"
#include "meta.h"

// dense containers
#include "vector.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "cmat3d.h"
#include "cmat4d.h"

// fixed size containers
#include "fixsize.h"

// sparse containers
#include "diag.h"
#include "matcoo.h"
#include "matcooh.h"
#include "cmatobd.h"

// dense slicing
#include "svector.h"
#include "scmat.h"

// stride slicing
#include "dvector.h"
#include "dcmat.h"
#include "jcmat.h"
#include "jcmat3d.h"
#include "jcmat4d.h"

// arithmetics
#include "imag.h"
#include "copy.h"
#include "complex_arith.h"
#include "scalar_arith.h"
#include "ptr_arith.h"
#include "arithmetic.h"
#include "arithmetic1.h"
#include "sparse_arith.h"
#include "slice_arith.h"

// special
#include "random.h"
#include "interp1.h"
#include "lin_eq.h"
#include "eig.h"
#ifdef SLS_USE_GSL
#include "ylm.h"
#endif
#include "anglib.h"
#include "coulomb.h"
#include "mat_fun.h"
#include "expokit.h"
#include "fedvr.h"
#include "fft.h"
#include "interv.h"
// #include "mparith.h"

// other utilities
#include "sort.h"
#include "search.h"
#include "input.h"
#include "file.h"
#include "matt.h"
#include "disp.h" // see also print.cpp
#include "time.h"

// string/unicode
#include "string.h"
#include "parser.h"
#include "unicode.h"
