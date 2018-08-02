=== Introduction ===
This project is my C++ utilities based on Numerical Recipes 3ed. It typedefs all of the numerical types it uses, and defines some simple template clases for vectors, matrices, and 3D matrices. All other utilities uses these types and classes.

"nr3plus.h" and "nr3plus.cpp" include the minimal set of this project, and can be used independently without other files.

Vectors or matrixs temporary is hard to deal with, so don't use it! This means, don't ever return a vector or matrix from a function (including operator overloading). Note that +=, -= etc. is totally fine because they don't return anything.


=== typedefs ===
Int is int, Uint is unsigned int, Llong is 64 bit integer, Doub is double, Char is char, Uchar is unsigned char, Ldoub is long double, 
Long is 64 bit integer by default, define _USE_Int_AS_LONG macro to use int instead.

Functions output arguments can not be any of the input arguments.

=== constants ===
Doub PI{ 3.14159265358979323 }
Complex I(0, 1)
Doub E{exp(1.)}

=== time utilities ===
void tic()
Doub toc()
void tic(Int ind)
Doub toc(Int ind)

=== scalar utilities ===
Int isodd(Int n) // return 1 if n is odd, return 0 otherwise
Bool ispow2(Int n) // if n is a power of 2 or 0
operator +,-,*,/ between Complex and Int

=== vec/mat display ===
void disp(av) // can also be used while debugging, because of this, default arguments are not allowed
void disp(v, start, n)
void disp(a, start1, start2, n1, n2)
void disp(a3, start1, start2, start3, n1, n2, n3)
void disp(..., precision)

=== get vec/mat properties ===
n = numel(av)
p = pointer(av) // get the pointer to the first element
s = sum(av)
s = max(av) // max of absolute values for complex mat/vec
s = max(ind, av) // also output the index
s = norm(v) // vec/mat norm
s = norm2(v) // vec/mat norm squared

=== matrix manipulation ===
a << s  // set vec/mat to a constant value
void linspace(vec/mat, first, last, N = -1)
void trans(a) // matrix transpose
void her(a) // hermitian conjugate
void flip(NRvector<T>)
void shift(NRmatrix<T> &a, const Int nshift, const Int dim = 1)
void diagonals(NRmatrix<T> &a) // shift the i-th line i times to the left, moving diagonals to columns
void idiagonals(NRmatrix<T> &a) // inverse of diagonals(), shift the i-th line i times to the right

=== vectorized math functions ===
sin(), cos(), exp(), tan(), whenever make sense

=== matrix arithmatics ===
operators +,-,*,/ scalar/vec/mat, whenever make sense (inefficient!).
operators +=,-=,*=,/= scalar/vec/mat, whenever make sense
void plus(out, in, in) for scalar/vec/mat, whenever make sense.
void minus(out, in, in)
void emul(out, in, in) // element-wise multiplication
void real(MatDoub_O &rc, MatComplex_I &c)
void imag(MatDoub_O &ic, MatComplex_I &c)
void abs(out, in), whenever make sense
void complex(VecComplex_O &y, VecDoub_I &x)
void conjugate(VecComplex_IO v)
s = v1*v2 // dot product, whenever make sense
void outprod(MatComplex_O &prod, VecComplex_I &v1, VecComplex_I &v2) // outer product
void mul(out, in, in) // mat-mat or mat-vec multiplications, whenever make sense

=== FT related ===
fftshift()
void dft(MatComplex_O &Y, Doub kmin, Doub kmax, Int Nk, MatComplex_I &X, Doub xmin, Doub xmax)
void idft(MatComplex_O &X, Doub xmin, Doub xmax, Int Nx, MatComplex_I &Y, Doub kmin, Doub kmax)
