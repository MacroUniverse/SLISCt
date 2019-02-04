Scientific Library In Simple C++ (SLISC)

## Introduction

SLISC is a header-only library mostly based on Numerical Recipes 3ed, using simple C++ 11 gramars so that it is easy to read, use and modify while maintaining a high performance. The library currencly provides simple class templates for vector, matrix and 3D matrix. Algorithms from Numerical Recipes are slowly added. The library also provides some utilities frequently used, such as timers and debug utilities. The library can optionally use some algorithms in Eigen (another header-only library).

A simple example :

```cpp
#include "SLISC/algorithm.h"
#include "SLISC/disp.h"
int main()
{
	using namespace slisc;
	using std::cout; using std::endl;
	VecDoub u(3), v(3); // vectors, double type
	linspace(u, 0, 2); // elements linearly spaced from 0 to 2
	cout << "u = \n"; disp(u); // display (print) vector/matrix
	v = 3.14; // set all elements to 3.14
	u += v; // vector plus vector
	v += 2; // vector plus scalar
	MatDoub a, b(1, 1); // matrices, double  type, always row major
	b.resize(2, 3); // resize b to 2 columns and 3 rows
	a.resize(b); // resize a to have the size of b
	a(0, 0) = 1.1; // access element by row and column indices
	a(3) = 9.9; // access element by a single index
	a.end() = 5.5; // last element
	cout << "a has " << a.nrows() << " rows and " << a.ncols()
	<< " columns, and a total of " << a.size() << " elements." << endl;
	disp(a);
}
```

SLISC has a modular design like the Standard Template Library. Just include any header file in SLISC folder. All definitions has namespace "slisc".

## programming style

* Class object temporary is inefficient (even with move constructor/assignment), using copy/move constructor or move assignment operator for vector/matrix types will create an error. Vector/Matrix type arguments should be passed by reference and should not be returned (use reference for output).

* Avoid using unsigned integer types as much as possible (this is also the google c++ style).

* Place a dot after Doub literals is prefered.

* Generally, functions output arguments can not be any of the input arguments (this is called aliasing).

* Intrinsic types are aliased inside the library. "Int" is 32-bit integer; "Uint" is "unsigned Int"; "Llong" is 64-bit integer; "Doub" is double (64-bit); "Comp" is "std::complex\<Doub>"; "Char" is "char"; "Uchar" is "unsigned char"; "Ldoub" is "long double"; "Long" is "Llong" by default, define "_USE_Int_AS_LONG" macro to use "Int" instead. Use "Long" for vector/matrix index or loop variable etc.

* A type ending with "_I" is the const (or reference to const) version of that type, used in function argument declarations to indicate input argument. Similarly, "_O" means output (reference type), "_IO" means both input and output (reference type), both of them are just the non-const version of the type. Note that a reference to "_O" or "_IO" types is still a reference type.

* Class members variables should start with `m_` for clearity, and avoid name confliction with member function arguments.

* Only use up to c++11 features.

* Will not use any serious meta-programming, to keep the code readable and easy to modify (even this means creating blocks of similar code).

## "slisc.h"
"slisc.h" includes some constants, type alias, and vector/matrix class template definitions.

```cpp
const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);
```

### Common Members for Vector, Matrix, Mat3d Templates

Operator () : get a reference for the i-th element, in row-major order.
Operator << : move data; empty vector/matrix can not be moved.

end() : get a reference for the last element.

size() : get the number of elements

### Vector Class Template

Constructors: Vector() for default, Vector(Long_I n) for vector size, Vector(Long_I n, const T &a) to specify element as well, Vector(Long_I n, const T *a) to initialize from array.

Operator = : Copy-assignment operator has auto resize, self-assignment is forbidden. The right hand side can be a scalar.

Operator [] : for vector, get a reference for the i-th element; for matrix, return a pointer for the second index.


resize(Long_I) : resize vector, contents are not preserved. resize() does nothing if size doesn't change.

resize(Vector<> v) : resize to the same size of v

### Matrix Class Template
The matrix template name is Matrix<T>, Methods are similar to that of vector class. Matrix is row-major only. 

TODO.

### Mat3d Class Template
3D matrix template name is Mat3d. Methods are similar to that of vector class. Mat3d is row-major only.

TODO.

### Vector/Matrix Type Alias
The typedefs for vector/matrix classes are (each type also comes with "_I", "_O", and "_IO" versions) :  VecInt, VecUint, VecLlong, VecUllong, VecChar, VecUchar, VecDoub, VecComp, VecBool, MatInt, MatUint, MatLlong, MatUllong, MatChar, MatUchar, MatDoub, MatComp, MatBool, Mat3Doub, Mat3Comp.

## algorithm.h
includes basic arithmatics like "==", "+=", "*=", plus(), minus(), etc.

## disp.h
includes various overloaded "disp()" functions.

## time.h
time utilities

// all times are in seconds.
void Timer::tic()
Doub Timer::toc()
void CPUTimer::tic()
Doub CPUTimer::toc()

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

If you want to use "disp()" in debugger, add "SLISC/print.cpp" to compiler and use "print()" with the same arguments.

## get vec/mat properties

```cpp
n = numel(av)
s = sum(av)
s = max(av) // max of absolute values for complex mat/vec
s = max(ind, av) // also output the index
s = norm(v) // vec/mat norm
s = norm2(v) // vec/mat norm squared
```

## matrix manipulation
```cpp
void linspace(vec/mat, first, last, N = -1)
void trans(a) // matrix transpose
void her(a) // hermitian conjugate
void flip(Vector<T>)
void shift(Matrix<T> &a, const Int nshift, const Int dim = 1)
void diagonals(Matrix<T> &a) // shift the i-th line i times to the left, moving diagonals to columns
void idiagonals(Matrix<T> &a) // inverse of diagonals(), shift the i-th line i times to the right
```

## element-wise math functions

sin(), cos(), exp(), tan(), whenever make sense


## matrix arithmatics

```cpp
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
```

## "print.cpp"
For debug purpose only. Need to compile "print.cpp" along with the program.

"print()" are same as "disp()" except they can be called in debugger. To realize this, print() does not have any namespace, and cannot be inlined.

In gdb, things are a little more complicated, since gdb does not fully support function overloading, use "print2()" for 2 arguments, "print3()" for 3 arguments, etc.

## Calculus

void integral(Vector<T> &F, const Vector<T> &f, Doub_I dx) // simple indefinite integration


## FFT related
```cpp
fftshift()
void dft(MatComplex_O &Y, Doub kmin, Doub kmax, Int Nk, MatComplex_I &X, Doub xmin, Doub xmax)
void idft(MatComplex_O &X, Doub xmin, Doub xmax, Int Nx, MatComplex_I &Y, Doub kmin, Doub kmax)
```

## String related

template\<typename T> inline std::string num2str(T s) // mainly std::to_string(), but no trailing zeros.

## OpenMP functions
```cpp
// parallelized version of functions
void diagonals_par(Matrix<T> &a)
void idiagonals_par(Matrix<T> &a)
void outprod_par(Matrix<T> &prod, const Vector<T1> &v1, const Vector<T2> &v2)
void outprod_par(Matrix<T> &prod, VecComp_I &v1, const Vector<T2> &v2)
void mul_par(Vector<T> &y, const Vector<T1> &x, const Matrix<T2> &a)
void dft_par(MatComp_O &Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I &X, Doub xmin, Doub xmax)
void idft_par(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax)
```

## Related Projects
* See cuSLISC project for a GPU version of SLISC using CUDA.
* See MatFile project for saving and reading "Vector" or "Matrix" to/from Matlab data file ".mat", or text based file ".matt".

## TODO
* I should define "I" as a spetial class, and implement more efficient "+", "-", "*", "/", etc.
* replace `error()` macro with `throw()`
* incorporate "arb" library for evaluation of some special functions, and for multi-precision arithmetic (does not work for windows yet)
* use BLAS/LAPACK to enhance performance (optionally), and time different implementations (mine, Eigen, MKL)
* implement column major matrix classes
* Add the "a(i,j)" format of matrix indexing for row majored matrix (and maybe consider abandoning "a[i][j]" format later).
* put all internal names into "slisc::internal" namespace
* test `randInt()`
* test single indexing operator() for matrix and mat3d
* update from Go-Solver project
* test McooH and arithmetic
* test "expokit.h"
* test operator= for different types
* test Comp2Real
* 