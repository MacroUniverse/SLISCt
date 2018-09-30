# Makefile

exe = nr3.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = \
main.cpp \
nr3plus.cpp \
fourier.cpp \
interp_1d.cpp \
interp_2d.cpp \
ludcmp.cpp \
tridag.cpp \
eigen_basics.cpp \
eigen_fft.cpp \
eigen_linsolve.cpp

objects = $(source:.cpp=.o)

compiler = g++

flags =  \
-I $(eigenPath)\
-std=c++11
# -g
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects) $(xobjects) $(libs)

$(objects):$(source)
	$(compiler) -c $(flags) $(source) $(sourcepath)$(xsource)

# must specify shared library directory
#run:
#	LD_LIBRARY_PATH=$(libpath) ./$(exe)

# clean all except source
clean:
	rm -f *.o *.x *.gch
