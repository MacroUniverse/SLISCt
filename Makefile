# Makefile

exe = slisc.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = main.cpp #SLISC/print.cpp
objects = main.o #print.o

compiler = icpc
# (NOTE:  The  icpc command uses the same compiler options as the icc com-
# mand. Invoking the compiler using icpc compiles .c and .i files as C++.
# Invoking  the  compiler  using icc compiles .c and .i files as C. Using
# icpc always links in  C++  libraries.  Using  icc  only  links  in  C++
# libraries if C++ source is provided on the command line.)

flags =  -I $(eigenPath) -std=c++17 -g -mkl -fp-model precise -fp-model except
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects)

$(objects):$(source)
	$(compiler) -c $(flags) $(source)

clean:
	rm -f *.o *.x *.gch
