# Makefile

exe = main.x

eigenPath = ../EigenTest/Eigen/

compiler = icpc

# use `sudo apt install libgsl-dev` to install GNU scientific library
# use `dpkg -L dpkg -L libgsl-dev` to check the installation directory
lib = -lgsl

flags = -Wall -I $(eigenPath) -std=c++14 -O3 -D NDEBUG -fopenmp $(no_warn) -mkl -fp-model precise -fp-model except -qopenmp $(lib)
# -g # debug
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
