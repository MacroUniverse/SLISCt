Scientific Library In Simple C++ (SLISC)

## Introduction

SLISC is a header-only library mostly based on Numerical Recipes 3ed, using simple C++ features so that it is easy to read, use and modify while maintaining a high performance. The library currencly provides simple class templates for vector, matrix (row-major and col-major, fixed-size and sparse) and 3D matrix. Algorithms from Numerical Recipes are slowly added. The library also provides some utilities frequently used, such as timers and IO utilities. The library can optionally use some algorithms in Eigen (another header-only library) or use Intel MKL subroutines.

A simple example :

```cpp
#include "SLISC/algorithm.h"
#include "SLISC/disp.h"
int main()
{
	using namespace slisc;
	VecDoub u(3), v(3); // vectors, double type
	linspace(u, 0, 2); // elements linearly spaced from 0 to 2
	cout << "u = \n"; disp(u); // display vector/matrix
	v = 3.14; // set all elements to 3.14
	u += v; // vector plus vector
	v += 2; // vector plus scalar
	MatDoub a, b(1, 1); // row major matrices of double precision
	b.resize(2, 3); // resize b to 2 columns and 3 rows
	a.resize(b); // resize a to have the size of b
	a(0, 0) = 1.1; // access element by row and column indices
	a[3] = 9.9; // access element by a single index
	a.end() = 5.5; // last element
	cout << "a has " << a.nrows() << " rows and " << a.ncols()
	<< " columns, and a total of " << a.size() << " elements." << endl;
	disp(a);
}
```

SLISC has a modular design like the Standard Template Library. Just include any header file in the `SLISC/` folder. All definitions has namespace `slisc`.

## Programming Style

* All containers types are returned by reference.
* Avoid using unsigned integer types when possible. They are not supported by SLISC for now.
* Generally, functions output arguments can not be any of the input arguments (this is called aliasing).

* Intrinsic types are aliased inside the library. For example, `Bool` is `bool`, `Int` is 32-bit integer, `Doub` is `double` (64-bit); `Comp` is `std::complex<Doub>`. `Long` is used as vector/matrix index and variable. `Long` is `Llong` by default, define `SLS_USE_INT_AS_LONG` macro to use "Int" instead. 

* A type with `_I` suffix is the `const` or `reference to const` version of that type, used in function parameter declarations to indicate an input argument. Similarly, `_O` means output (reference type), `_IO` means both input and output (reference type). Note that a reference to `_O` or `_I` types is still a reference type.

## Headers Introduction
* `slisc.h` includes all type definitions and constants, and basic arithmetics.
* `global.h` has all the container declaration and type definitions etc.
* `meta.h` has all the meta-programming utilities.
* `complex_arith.h` defines extra operators involving std::complex<>, such as  `+, -, *, /, +=, -=, *=, /=, ==, !=`.
* `scalar_arith.h` defines scalar utilities such as `MIN()`, `MAX()`, `SQR()`, `isodd()`, `mod()`.
* `vector.h` defines the base type `Vbase<T>` and vector container `Vector<T>`.
* `matrix.h` defines the row-major matrix container `Matrix<T>`.
* `cmat.h` defines the col-major matrix container `Cmat<T>`.
* `fixsize.h` defines the fixed-size vector `FixVec<T,N>`, and col-major matrix `FixCmat<T,Nr,Nc>`.
* `mat3d.h` defines the row-major 3D array `Mat3d<T>`.
* `sparse.h` defines the sparse square diagonal matrix `Diag<T>`, COO sparse matrix `MatCoo<T>`, COO sparse Hermitian matrix `MatCooH<T>`.
* TODO...

## Scalar Types
TODO...
For example, `Bool` is `bool`, `Char` is `char`, `Int` is 32-bit integer, `Llong` is 64-bit integer; `Float` is `float`, `Doub` is `double` (64-bit); `Comp` is `std::complex<Doub>`; `Ldoub` is `long double`; `Long` is `Llong` by default, define "SLS_USE_INT_AS_LONG" macro to use "Int" instead. Use "Long" for vector/matrix index or loop variable etc.

## Constants
```cpp
const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);
```

### Common Members for Vector, Matrix, Mat3d Templates
* `size()`: return total number of elements.
* `nrows(), ncols()`: return number of rows or columns.
* `operator[][i], operator()(i)` : return a reference for the i-th element.
* `operator()(i,j), (i,j,k)` : return a reference for an element.
* `end()` : return a reference for the last element.
* `operator<<` : transfer data to another container.
* TODO...

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
* includes basic arithmatics like "==", "+=", "*=", plus(), minus(), etc.
* Operators `+, -, *, /, +=, -=, *=, /=` are only for container types element-wise operations.

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

* `sin(v, v1)`, `cos(v, v1)`, `exp(v, v1)`, `tan(v, v1)`, etc.

## matrix arithmatics

```cpp
operator ==, != // compare size and each element, right hand side can also be a scalar.
operators +,-,*,/ scalar/vec/mat, whenever make sense (inefficient!).
operators +=,-=,*=,/= scalar/vec/mat, whenever make sense
void Plus(out, in, in) //for scalar/vec/mat, whenever make sense.
void Minus(out, in, in) // binary "-" operator
void Minus(in_out) // unary "-" operator
void Minus(out, in)
Times, Divide
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

```cpp
void integral(Vector<T> &F, const Vector<T> &f, Doub_I dx) // simple indefinite integration
```

## FFT related
```cpp
fftshift()
void dft(MatComplex_O &Y, Doub kmin, Doub kmax, Int Nk, MatComplex_I &X, Doub xmin, Doub xmax)
void idft(MatComplex_O &X, Doub xmin, Doub xmax, Int Nx, MatComplex_I &Y, Doub kmin, Doub kmax)
```

## String related
```cpp
template<typename T> inline Str num2str(T s) // mainly std::to_string(), but no trailing zeros.
```

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

## Data File (.matt)
.matt is a text based format that I designed to immitate Matlab's .mat file.
* can write many variables and matrices to a single file, with names.
* can find and read any number of stored variables/matrices from a file by name without scanning through the whole file. Thus, reading a large file is fast.
* there is a Matlab function `mattread()` that works similar to `load()` to read this file into Matlab workspace. This reading subroutine can be easily implemented for any other language.
* Matrices are stored in column major order. Memory access speed is negaligible comparing to hard-drive.
* TODO: is it possible to use space instead of `\n` so that there is no CRLF issue?
* TODO: use `typenum<T>()` to specify type, instead of 1,2,3,etc.
* TODO: `i` is not necessary, consider use `a+b` or `a-b` to represents complex numbers. However, `mattread()` might take longer time because Matlab cannot recognize this format directly.

## Related Projects
* See cuSLISC project for a GPU version of SLISC using CUDA.
* See MatFile project for saving and reading "Vector" or "Matrix" to/from Matlab data file ".mat", or text based file ".matt".

## Internal Coding Rules
* C++14/17 features used: `if constexpr`
* Code should be easy to understand and modify.
* Interface and implementation should separated.
* Class members variables should start with `m_` for clearity, and avoid name confliction with member function arguments.
* Use SFINAE macro `SLISC_IF(bool)` to limit template instanciation.
* Templates must work for all possible instanciations.

## TODO
* consider define pure imaginary number as a spetial class, and implement more efficient "+", "-", "*", "/", etc.
* consider replacing `error()` macro with `throw()`
* `warning()` macro is already defined!
* incorporate "arb" library for evaluation of some special functions, and for multi-precision arithmetic (does not work for windows yet)
* use MKL to enhance performance (optionally), and time different implementations (mine, Eigen, MKL)
* put all internal names into "slisc::internal" namespace
* test `randInt()`
* update from Go-Solver project
* test "expokit.h"
* test sparse arithmetic
* ptr_arith.h functions must support N = 0
* complete test of containers
* complete test of "meta.h"
* complete test of "ptr_arith.h"
* complete test of "arithmetic.h"
* check all compiler warnings
* use `SLS_FOTBID_COPY_CONSTRUCTOR` to forbit copy constructor of containers, default should be undefined.
* test "meta.h" for ImagNum<T> types
