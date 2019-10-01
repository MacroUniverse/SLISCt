# choose one of the following makefiles for different compilers and libraries
# edit macros in global.h in correspondence

# macros: SLS_USE_CBLAS SLS_USE_LAPACKE
# known bug: LAPACKE is not working for now
include make/compatible.mak

# macros: SLS_USE_MKL 
# include make/g++_and_mkl.mak

# macros: SLS_USE_MKL
# include make/icpc_and_mkl.mak
