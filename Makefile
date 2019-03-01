# Makefile

exe = slisc.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = main.cpp #SLISC/print.cpp
objects = main.o #print.o

compiler = icc

flags =  -I $(eigenPath) -std=c++17 -g -mkl -fp-model precise -fp-model except
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects)

$(objects):$(source)
	$(compiler) -c $(flags) $(source)

clean:
	rm -f *.o *.x *.gch
