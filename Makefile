# Makefile

exe = slisc.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = main.cpp SLISC/disp.cpp
objects = main.o disp.o

compiler = g++

flags =  \
-I $(eigenPath) \
-std=c++11 \
-g
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects) $(xobjects) $(libs)

$(objects):$(source)
	$(compiler) -c $(flags) $(source) $(sourcepath)$(xsource)

clean:
	rm -f *.o *.x *.gch
