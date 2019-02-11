// all containers and basic utilities
// also list all dependencies in order
#pragma once

// === user macros ===
// #define SLS_CHECK_BOUNDS // defined when debugging
// #define SLS_CHECK_SHAPE  // defined when debugging

// #define SLS_RAND_SEED 1 // defined as std::time()
// #define SLS_USE_INT_AS_LONG // undefined
// #define SLS_CHECK_COO_REPEAT // defined when debugging

// === internal macros ===
// error(str)
// warning(str)
// SLS_IF_HELPER(cond)
// SLS_IF(...)
// time_h_error(str)

#include "global.h"
#include "meta.h"
#include "complex_arith.h"
#include "scalar_arith.h"

// containers
#include "vector.h"
#include "matrix.h"
#include "cmat.h"
#include "mat3d.h"
#include "fixsize.h"
#include "ptr_arith.h"
#include "arithmetic.h"

#include "sparse.h"
#include "sparse_arith.h"

#include "disp.h"
