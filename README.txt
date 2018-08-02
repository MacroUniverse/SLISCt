=== Introduction ===
This project is a simple C++ math library rewritten from Numerical Recipes 3ed. It typedefs some common scalar types, and some simple template clases for vectors, matrices, and 3D matrices. All other utilities uses these types and classes.

"nr3.h" includes all the typedefs and vector/matrix class definitions and dependencies, and can be used independently;

"nr3plus.h" and "nr3plus.cpp" include the common utilities, and only depends on "nr3.h".

Class object temporary is inefficient (even with move constructor/assignment), using copy/move constructor or move assignment operator for vector/matrix types will create an error. Vector/Matrix type arguments should be passed by reference and should not be returned (use reference for output).

=== typedefs ===
Int is int, Uint is unsigned int, Llong is 64 bit integer, Doub is double, Comp is std::complex<double>, Char is char, Uchar is unsigned char, Ldoub is long double.
Long is 64 bit integer by default, define _USE_Int_AS_LONG macro to use int instead. Use Long for vector/matrix index or loop variable etc.
Types ending with "_I" is const version of that type, used in function argument declarations to indicate input argument. Similarly, "_O" means output, "_IO" means both input and output, both of them are just the non-const version of the type.

By default, functions output arguments can not be any of the input arguments (this is called aliasing).

=== constants ===
const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);

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


=== others ===
resize() does nothing if size doesn't change
