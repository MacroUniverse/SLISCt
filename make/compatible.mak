# use g++ with libgsl-dev libblas-dev liblapacke-dev (all available from apt-get)

# use `dpkg -L dpkg -L libgsl-dev` to check the installation directory
libs = -lgsl -llapacke -lblas

flags = -Wall -Wno-reorder -fopenmp

# link
# choose `$(mkl_dyn_link)` or `$(mkl_stat_link)`
goal:main.o
	g++ $(flags) -o main.x main.o $(libs)

# compile
main.o:main.cpp
	g++ $(flags) -std=c++17 -c main.cpp

clean:
	rm -f *.o *.x
