Scientific Library In Simple C++ (SLISC)

## Introduction

This project is a scientific library rewritten from Numerical Recipes 3ed, using simple C++ gramars so that it is easy to read, use and modify, while maintaining a high performance. The library currencly provides class templates for vector, matrix and 3D matrix, and some types that has fixed size in memory. Basic matrix/vector manipulation is provided. The library also has some utilities frequently used, such as timers, debug utilities. The library uses standard C++11, and has only 3 files, no other dependency is needed.

A simple example :

```cpp
#include "slisc.h" // this is now named "nr3plus.h"
using std::cout; using std::endl;
int main()
{
	VecDoub u(3), v(3); // vectors, double type
	linspace(u, 0, 2); // elements linearly spaced from 0 to 2
	cout << "u = \n"; disp(u); // display (print) vector/matrix
	v = 3.14; // set all elements to 3.14
	u += v; // vector plus vector
	v += 2; // vector plus scalar
	MatDoub a, b(1, 1); // matrices, double  type, always row major
	b.resize(2, 3); // resize b to 2 columns and 3 rows
	a.resize(b); // resize a to have the size of b
	a[0][0] = 1.1; // access element by row and column indices
	a(3) = 9.9; // access element by a single index
	a.end() = 5.5; // last element
	cout << "a has " << a.nrows() << " rows and " << a.ncols()
	<< " columns, and a total of " << a.size() << " elements." << endl;
	disp(a);
}
```

"nr3.h" includes all the typedefs and vector/matrix class definitions and dependencies, and can be used independently.

"nr3plus.h" and "nr3plus.cpp" include the rest of the library, and only depends on "nr3.h".


## programming style

Class object temporary is inefficient (even with move constructor/assignment), using copy/move constructor or move assignment operator for vector/matrix types will create an error. Vector/Matrix type arguments should be passed by reference and should not be returned (use reference for output).

Avoid using unsigned integer types as much as possible (this is also the google c++ style).

Place a dot after Doub literals is prefered.

=== typedefs ===
Int is int, Uint is unsigned int, Llong is 64-bit int, Doub is double, Comp is std::complex<double>, Char is char, Uchar is unsigned char, Ldoub is long double.
Long is 64 bit integer by default, define _USE_Int_AS_LONG macro to use int instead. Use Long for vector/matrix index or loop variable etc.
Types ending with "_I" is const version of that type, used in function argument declarations to indicate input argument. Similarly, "_O" means output, "_IO" means both input and output, both of them are just the non-const version of the type.

Generally, functions output arguments can not be any of the input arguments (this is called aliasing), except for element-wise functions.


## vector/matrix class template

Constructors: NRvector() for default, NRvector(Long_I n) for vector size, NRvector(Long_I n, const T &a) to specify element as well, NRvector(Long_I n, const T *a) to initialize from array.
Operator = : Copy-assignment operator has auto resize, self-assignment is forbidden. The right hand side can be a scalar.
Operator [] : for vector, get a reference for the i-th element; for matrix, return a pointer for the second index.
Operator () : get a reference for the i-th element, in row-major order.
Operator << : move data; empty vector/matrix can not be moved.
end() : get a reference for the last element.
size() : get the number of elements
resize(Long_I) : resize vector, contents are not preserved. resize() does nothing if size doesn't change.
resize(NRvector<> v) : resize to the same size of v

The matrix template name is NRmatrix<T>, 3D matrix template name is NRMat3d. Matrix is row-major only. Methods are similar to that of vector class.

The typedefs for vector/matrix classes are (each type also comes with "_I", "_O", and "_IO" versions) :  VecInt, VecUint, VecLlong, VecUllong, VecChar, VecUchar, VecDoub, VecComp, VecBool, MatInt, MatUint, MatLlong, MatUllong, MatChar, MatUchar, MatDoub, MatComp, MatBool, Mat3Doub, Mat3Comp.

NRbase class should never be used.

## constants

const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);


## time utilities

// all times are in seconds.
void tic()
Doub toc()
void tic(Int ind)
Doub toc(Int ind)
void ctic() // cpu time
Doub ctoc()


## scalar utilities

Int isodd(Int n) // return 1 if n is odd, return 0 otherwise
Bool ispow2(Int n) // if n is a power of 2 or 0
operator +,-,*,/ between Complex and Int


## vec/mat display

void disp(av) // can also be used while debugging, because of this, default arguments are not allowed
void disp(v, start, n)
void disp(a, start1, start2, n1, n2)
void disp(a3, start1, start2, start3, n1, n2, n3)
void disp(..., precision)


## get vec/mat properties

n = numel(av)
p = pointer(av) // get the pointer to the first element
s = sum(av)
s = max(av) // max of absolute values for complex mat/vec
s = max(ind, av) // also output the index
s = norm(v) // vec/mat norm
s = norm2(v) // vec/mat norm squared


## matrix manipulation

void linspace(vec/mat, first, last, N = -1)
void trans(a) // matrix transpose
void her(a) // hermitian conjugate
void flip(NRvector<T>)
void shift(NRmatrix<T> &a, const Int nshift, const Int dim = 1)
void diagonals(NRmatrix<T> &a) // shift the i-th line i times to the left, moving diagonals to columns
void idiagonals(NRmatrix<T> &a) // inverse of diagonals(), shift the i-th line i times to the right


## element-wise math functions

sin(), cos(), exp(), tan(), whenever make sense


## matrix arithmatics

operator ==, != // compare size and each element, right hand side can also be a scalar.
operators +,-,*,/ scalar/vec/mat, whenever make sense (inefficient!).
operators +=,-=,*=,/= scalar/vec/mat, whenever make sense
void plus(out, in, in) //for scalar/vec/mat, whenever make sense.
void minus(out, in, in) // binary "-" operator
void minus(in_out) // unary "-" operator
void minus(out, in)
void emul(out, in, in) // element-wise multiplication
void real(MatDoub_O &rc, MatComplex_I &c)
void imag(MatDoub_O &ic, MatComplex_I &c)
void abs(out, in), whenever make sense
void complex(VecComplex_O &y, VecDoub_I &x)
void conjugate(VecComplex_IO v)
s = v1*v2 // dot product, whenever make sense
void outprod(MatComplex_O &prod, VecComplex_I &v1, VecComplex_I &v2) // outer product
void mul(out, in, in) // mat-mat or mat-vec multiplications, whenever make sense


## calculus

void integral(NRvector<T> &F, const NRvector<T> &f, Doub_I dx) // simple indefinite integration


## FT related

fftshift()
void dft(MatComplex_O &Y, Doub kmin, Doub kmax, Int Nk, MatComplex_I &X, Doub xmin, Doub xmax)
void idft(MatComplex_O &X, Doub xmin, Doub xmax, Int Nx, MatComplex_I &Y, Doub kmin, Doub kmax)


## string related

template <typename T> inline std::string num2str(T s) // mainly std::to_string(), but no trailing zeros.


## OpenMP functions

// parallelized version of functions
void diagonals_par(NRmatrix<T> &a)
void idiagonals_par(NRmatrix<T> &a)
void outprod_par(NRmatrix<T> &prod, const NRvector<T1> &v1, const NRvector<T2> &v2)
void outprod_par(NRmatrix<T> &prod, VecComp_I &v1, const NRvector<T2> &v2)
void mul_par(NRvector<T> &y, const NRvector<T1> &x, const NRmatrix<T2> &a)
void dft_par(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax)
void idft_par(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax)
