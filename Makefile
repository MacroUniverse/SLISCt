# Makefile

exe = slisc.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = main.cpp SLISC/print.cpp
objects = main.o print.o

compiler = g++

flags =  -I $(eigenPath) -std=c++11 -g
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects)

$(objects):$(source)
	$(compiler) -c $(flags) $(source)

clean:
	rm -f *.o *.x *.gch
