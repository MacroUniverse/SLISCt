#pragma once
#include "interp_1d.h"

// my function for interp. a complex grid from a complex grid
// Z[i][j] = z(y[i], x[j])
// x, y, x0, y0 must be in ascending order
// Z[i][j] = 0 for extrapolation
void Spline_grid(MatComp_O &Z, VecDoub_I &x, VecDoub_I &y, MatComp_I &Z0, VecDoub_I &x0, VecDoub_I &y0);

