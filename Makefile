# Makefile

exe = main.x

eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

compiler = icpc
# (NOTE:  The  icpc command uses the same compiler options as the icc com-
# mand. Invoking the compiler using icpc compiles .c and .i files as C++.
# Invoking  the  compiler  using icc compiles .c and .i files as C. Using
# icpc always links in  C++  libraries.  Using  icc  only  links  in  C++
# libraries if C++ source is provided on the command line.)

# use `sudo apt install libgsl-dev` to install GNU scientific library
# use `dpkg -L dpkg -L libgsl-dev` to check the installation directory
lib = -lgsl

flags = -Wall -I $(eigenPath) -std=c++17 -g -fopenmp $(no_warn) -mkl -fp-model precise -fp-model except -qopenmp $(lib)
# -O3 # highest optimization
# -qopenmp # run OpenMP in parallel mode
# -qopenmp-stubs # run OpenMP in serial mode
# -fp-model # floating point model

# goal
goal: main.o
	$(compiler) -o $(exe) $(flags) main.o

main.o: main.cpp SLISC/*.h test/*.h
	$(compiler) $(flags) -c main.cpp

clean:
	rm -f *.o *.x *.gch
