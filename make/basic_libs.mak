# use g++ with libgsl-dev libblas-dev liblapacke-dev (all available from apt-get)
# use `dpkg -L dpkg -L lib***` to check the installation directory

libs = -lgsl -llapacke -lblas

flags = -Wall -Wno-reorder -fopenmp -O3 -D NDEBUG -D SLS_USE_CBLAS -D SLS_USE_LAPACKE -D SLS_USE_GSL

compiler = g++

# link
# choose `$(mkl_dyn_link)` or `$(mkl_stat_link)`
goal:main.o
	$(compiler) $(flags) -o main.x main.o $(libs)

# compile
main.o:main.cpp
	$(compiler) $(flags) -std=c++14 -c main.cpp

clean:
	rm -f *.o *.x
