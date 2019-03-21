Scientific Library In Simple C++ (SLISC)

## Introduction

SLISC is a header-only library written in a style similar to Numerical Recipes 3ed, using simple C++ features so that it is easy to read and modify while maintaining a high performance. The library currencly provides simple class templates for vector, matrix (row-major and col-major, fixed-size and sparse), 3D matrix (row-major), and basic arithmetics for them. Codes from many other projects or libraries has been incorporated into SLISC (e.g. Numerical Recipes, Eigen, Intel MKL etc.). The library also provides some utilities frequently used, such as timers and IO utilities (a text based file format `.matt` similar to Matlab's `.mat`).

SLISC has a comprehensive test suit, main.cpp will execute all the tests. Tests has been performed in Windows using Visual C++ and Intel compilers in Visual Studio 15.9.8 and in Linux using gcc and Intel compilers. If intel MKL (now free) is not installed, some functions will not work or will be much slower. Note that SLISC is a project in development, interfaces are subjected to change and not all code are working.

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

SLISC has a modular design like the Standard Template Library. Just include any header file(s) in the `SLISC/` folder. All definitions has namespace `slisc`.

## Recommended Programming Style
* All SLISC containers types (e.g. Matrix<>, Vector<>) should be returned by reference.
* Avoid using unsigned integer types when possible. They are not supported by SLISC for now.
* Generally, functions output arguments can not be any of the input arguments (this is called aliasing).

* Intrinsic types are aliased inside the library. For example, `Bool` is `bool`, `Int` is 32-bit integer, `Doub` is `double` (64-bit); `Comp` is `std::complex<Doub>`. `Llong` is `long long`.

* A type with `_I` suffix is the `const` or `reference to const` version of that type, used in function parameter declarations to indicate an input argument. Similarly, `_O` means output (reference type), `_IO` means both input and output (reference type). Note that a reference to `_O` or `_I` types is still a reference type.

## Meta Programming
* The file "meta.h" implements some meta programming utilities. These utilities are as user-friendly as possible, so that anyone with a basic c++ knowledge can understand their meaning with a glance. SLISC users should have a basic idea of how these utilities work, and are welcome to use them, but there is no need to understand how they are implemented.

* The macro function `SLS_IF(condition)` can be used to limit template function instantiation. If `condition` is `true`, then template can be instantiated normally, if `condition` is `false` the function will not be instantiated. As an example, the function `mul()` in "arithmetic.h" for matrix-vector multiplication is defined as
```cpp
template <class T, class T1, class T2, SLS_IF(is_Vector<T>() && is_dense_mat<T1>() && is_Vector<T2>())>
inline void mul(T &y, const T1 &a, const T2 &x)
{/*...*/}
```
Where `is_Vector<T>()` returns `true` if `T` is a `Vector<>` container, `is_dense_mat<T>()` returns `true` if `T` is a dense matrix (such as `Matrix<>`, `Cmat<>`, `FixCmat<>`). These function templates are also defined in "meta.h".

## Headers Introduction
When using something in any header file, just including that header file will be enough. Header files can be included in any order. Here is some brief introduction for each header file:
* `slisc.h` is a shorthand for some of the mostly used headerfiles.
* `global.h` has all the container declaration and type definitions etc.
* `meta.h` has all the meta-programming utilities.
* `complex_arith.h` defines extra operators used by std::complex<>, such as  `+, -, *, /, +=, -=, *=, /=, ==, !=`.
* `imag.h` defines a pure imaginary number type `Imag<>`, and related operators.
* `scalar_arith.h` defines scalar utilities such as `MIN()`, `MAX()`, `SQR()`, `isodd()`, `mod()`.
* `vector.h` defines the base type `Vbase<T>` and vector container `Vector<T>`.
* `matrix.h` defines the row-major matrix container `Matrix<T>`.
* `cmat.h` defines the col-major matrix container `Cmat<T>`.
* `fixsize.h` defines the fixed-size vector `FixVec<T,N>`, and col-major matrix `FixCmat<T,Nr,Nc>`.
* `mat3d.h` defines the row-major 3D array `Mat3d<T>`.
* `disp.h` display SLISC containers (matrix, vector, etc.)
* `input.h` promp for input, can save input history and repeat input automatically
* `matt.h` save/load text-based data files in `.matt` format, can save multiple named scalars and containers to a single file
* `ptr_arith.h` low level functions for `arithmetic.h`, using pointers as input and output instead of vector/matrix containers.
* `arithmetic.h` has utilities for dense matrices and vectors, e.g. `sum()`, `norm()`, dot product, matrix-vector multiplication.
* `slice.h` (experimental) matrix slicing, e.g. separate one column of a matrix and name it as a vector.
* `random.h` random number utilities
* `time.h` timing utilities
* `sort.h` sorting utilities
* `search.h` search elements in containers
* `string.h` string utilities
* `svd.h` for singlar value decomposition
* `eig.h` calculate matrix eigen values/vectors
* `fft.h` for fourier transforms
* `interp1.h` for 1 dimensional interpolation
* `interp2.h` for 2 dimensional interpolation
* `ludcmp.h` for LU decomposition
* `sparse.h` defines the sparse square diagonal matrix `Diag<T>`, COO sparse matrix `MatCoo<T>`, COO sparse Hermitian matrix `MatCooH<T>`, and basic arithmetics.
* `mat_fun.h` functions of square matrix
* `anglib.h` has functions for Clebsch–Gordan coefficients, 3j, 6j, and 9j symbols.
* `coulomb.h` calculates coulomb functions (F, G, H), and their derivatives.
* `fedvr.h` utilities for Finite Element Discrete Variable Representations, could be used to solve TDSE.
* `flm.h` a data structure for quantum mechanics wave functions in spherical coordinates (partial waves)
* `mparith.h` for arbitrary precision calculation
* `expokit.h` calculate matrix exponential

* TODO...

## "global.h"

### Scalar Types
`Bool` is `bool`, `Char` is `char`, `Int` is 32-bit integer, `Llong` is 64-bit integer; `Float` is `float`, `Doub` is `double` (64-bit); `Comp` is `std::complex<Doub>`; `Ldoub` is `long double`; `Long` `Long` is used as vector/matrix indices variables, and is `Llong` by default, define `SLS_USE_INT_AS_LONG` macro to use `Int` as `Long`.
TODO...

### Constants
```cpp
const Doub PI = 3.14159265358979323;
const Doub E  = 2.71828182845904524;
const Comp I(0., 1.);
```

## "meta.h":
* Every supported scalar type has a type number. `type_num<T>()` will return the type number of type `T`.
* Functions like `is_*(...)` can dynamically or statically check properties of types. The inputs are type numbers.
* functions like `is_*<...>()` can statically check properties of types. template parameters are types.
* `is_same<T1,T2>()` checks if `T1 = T2`.
* `is_integral<T>()` checks if `T` is an integral type, i.e. character types or integer types.
* `is_arithmetic<T>()` checks if `T` is an arithmetic type, i.e. integral type or floating point types.
* `is_fundamental<T>()` checks if `T` is a fundamental type, i.e. arithmetic type or `void` or `nullptr`.
* `is_signed<T>` checks if a type is signed, i.e. `T(-1) < T(0)`.
* the following utilities checks if `T` is a specific scalar type: `is_Bool<T>()`, `is_Char<T>()`, `is_Uchar<T>()`, `is_Int<T>()`, `is_Llong<T>()`, `is_Float<T>()`, `is_Ldoub<T>()`, `is_Fcomp<T>()`, `is_Comp<T>()`, `is_Lcomp()`, `is_Imag()`.
* `is_real<T>()` checks if `T` is not a `std::complex<>` or `Imag<>` type.
* `is_comp<T>()` checks if `T` is an `std::complex<>` type.
* `is_imag<T>()` checks if `T` is an `Imag<>` type.
* `is_scalar<T>()` checks if `T` is scalar type supported by SLISC.	
* `is_Vector<T>()` checks if `T` is a `Vector<>` type. The following are similar: `is_Matrix<T>()`, `is_Cmat()`, `is_FixVec()`, `is_FixCmat()`, `is_Mat3d()`, `is_Diag()`, `is_MatCoo()`, `is_MatCooH()`.
* `is_fixed<T>()` checks if `T` is a fixed-size container.
* `is_dense_mat<T>()` checks if is dense matrix (2D).
* `is_dense<T>()` checks if is dense container (including fixed-size)
* `is_sparse<T>()` check if is sparse vector/matrix
* Every SLISC matrix/vector container should has a container number, which can be returned by `contain_num<T>()`.
* `contain_type<T>()` returns the value type of any SLISC container type `T`.
* `is_real_dense<T>()` checks if `is_dense<T>() && is_real<contain_type<T>>()`.
* `is_comp_dense<T>()` checks if `is_dense<T>() && is_comp<contain_type<T>>()`
* `is_real_contain<T>()` checks if `is_contain<T>() && is_real<contain_type<T>>()`
* `is_comp_contain<T>()` checks if `is_contain<T>() && is_comp<contain_type<T>>()`
* `is_same_contain<T>()` checks if `is_contain<T1>() && contain_num<T1>() == contain_num<T2>()`

TODO...

## "vector.h"

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

`operator=` : Copy-assignment operator has auto resize, self-assignment is forbidden. The right hand side can be a scalar.

`operator[]` : for vector, get a reference for the i-th element; for matrix, return a pointer for the second index.


`resize(Long_I)` : resize vector, contents are not preserved. resize() does nothing if size doesn't change.

`resize(Vector<> v)` : resize to the same size of v

## "matrix.h"
The matrix template name is `Matrix<>`, Methods are similar to that of vector class. Matrix is row-major only. 

TODO.

## "mat3d.h"
3D matrix template name is `Mat3d<>`. Methods are similar to that of vector class. Mat3d is row-major only.

TODO.

## "slice.h"
`Svector<>` inherits `Vector<>` and thus can be casted to a vector when input to a function. `Svector` does not have it's own allocated memory, but use a block of contiguous memory from other dense containers (this is called slicing). For example, if we need to calculate the sum of a column of a `Cmat`, we can create an `Svector` to represent one column of `Cmat`, then use it as input to `sum()` function (need to cast to `Vector<>` first, unless `sum()` accepts `Svector` directly).

### Vector/Matrix Type Alias
The typedefs for vector/matrix classes are (each type also comes with "_I", "_O", and "_IO" versions) :  VecInt, VecUint, VecLlong, VecUllong, VecChar, VecUchar, VecDoub, VecComp, VecBool, MatInt, MatUint, MatLlong, MatUllong, MatChar, MatUchar, MatDoub, MatComp, MatBool, Mat3Doub, Mat3Comp.

## algorithm.h
* includes basic arithmatics like "==", "+=", "*=", plus(), minus(), etc.
* Operators `+, -, *, /, +=, -=, *=, /=` are only for container types element-wise operations.

## "disp.h"
includes various overloaded "disp()" functions.

## "time.h"
time utilities

// all times are in seconds.
void Timer::tic()
Doub Timer::toc()
void CPUTimer::tic()
Doub CPUTimer::toc()

## "scalar_arith.h"

Int isodd(Int n) // return 1 if n is odd, return 0 otherwise
Bool ispow2(Int n) // if n is a power of 2 or 0
operator +,-,*,/ between Complex and Int

## "disp.h"

void disp(av) // can also be used while debugging, because of this, default arguments are not allowed
void disp(v, start, n)
void disp(a, start1, start2, n1, n2)
void disp(a3, start1, start2, start3, n1, n2, n3)
void disp(..., precision)

If you want to use "disp()" in debugger, add "SLISC/print.cpp" to compiler and use "print()" with the same arguments.

## "arithmetic.h"

### basic utilities
```cpp
s = sum(av) // sum of elements
s = max(av) // max of absolute values for complex mat/vec
s = max(ind, av) // also output the index
s = norm(v) // vec/mat norm
s = norm2(v) // vec/mat norm squared
```

### matrix manipulation
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

## "fft.h"
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
parallelized version of a function should always end with "_par"
```cpp
// arithmetic.h
template <class T, class T1, class T2, SLS_IF(is_dense_mat<T>() && is_Vector<T1>() && is_Vector<T2>())>
inline void outprod_par(T &v, const T1 &v1, const T2 &v2);

template <class T, class T1, class T2, SLS_IF(is_Vector<T>() && is_dense_mat<T1>() && is_Vector<T2>())>
inline void mul_par(T &y, const T1 &a, const T2 &x);

template <class T, class T1, class T2, SLS_IF(is_Vector<T>() && is_Vector<T1>() && is_dense_mat<T2>())>
inline void mul_par(T &y, const T1 &x, const T2 &a)

// arithmetic1.h
void diagonals_par(Matrix<T> &a);
void idiagonals_par(Matrix<T> &a);

// fft.h
void dft_par(MatComp_O Y, Doub kmin, Doub kmax, Long_I Nk, MatComp_I X, Doub xmin, Doub xmax);
void idft_par(MatComp_O &X, Doub xmin, Doub xmax, Long_I Nx, MatComp_I &Y, Doub kmin, Doub kmax);
```

## Matt File (.matt)
* .matt is a text based format that I designed to immitate Matlab's .mat file.
* can write many variables and matrices to a single file, with names.
* can find and read any number of stored variables/matrices from a file by name without scanning through the whole file. Thus, reading a large file is fast.
* the read in value type T2 can be different from the written value type T1, if T1 can be losslessly converted to T2 (see `is_promo<T1,T2>()`).
* there is a Matlab function `save()` that works similar to `load()` to read this file into Matlab workspace. This reading subroutine can be easily implemented for any other language.
* Matrices are stored in column major order. Memory access speed is negaligible comparing to hard-drive.

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
* a constructor of Vbase/Vector that leaves things uninitialized might be added, and use it to optimize Svector() constructor
* modify "meta.h" so that `Svector` could be used as function arguments without casting to `Vector` first.
* modify `resize()` of dense matrices so that the input can be a sparse matrix
