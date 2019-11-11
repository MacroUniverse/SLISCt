# use g++ with MKL

# reference: MKL link advisor https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

mkl_stat_link = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

mkl_dyn_link = -L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

mkl_compile =  -m64 -I${MKLROOT}/include

# use `sudo apt install libgsl-dev` to install GNU scientific library
# use `dpkg -L dpkg -L libgsl-dev` to check the installation directory
libs = -lgsl

no_warn = -Wno-reorder # supress warning about initializer list order

flags = -Wall $(no_warn) -fopenmp -O3 -D NDEBUG -D SLS_USE_GSL -D SLS_USE_MKL

# link
# choose `$(mkl_dyn_link)` or `$(mkl_stat_link)`
goal:main.o
	g++ $(flags) -o main.x main.o $(libs) $(mkl_stat_link)

# compile
main.o:main.cpp
	g++ $(flags) $(mkl_compile) -std=c++14 -c main.cpp

clean:
	rm -f *.o *.x
