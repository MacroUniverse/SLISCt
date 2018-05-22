#pragma once
#include "interp_1d.h"

// my function for interp. a complex grid from a complex grid
// z[i][j] = z(y[end-i], x[j])
// x, y, x0, y0 must be in ascending order
// Z[i][j] = 0 for extrapolation
void Spline_grid(MatComplex_O &Z, VecDoub_I &x, VecDoub_I &y, MatComplex_I &Z0, VecDoub_I &x0, VecDoub_I &y0);

