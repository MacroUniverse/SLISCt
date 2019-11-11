# choose one of the following makefiles for different compilers and libraries

# the most compatible build, no libraries are required
# include make/compatible.mak

# blas lapacke gsl are required
include make/basic_libs.mak

# MKL and gsl are required
# include make/g++_and_mkl.mak

# icpc compiler and MKL are required
# include make/icpc_and_mkl.mak
