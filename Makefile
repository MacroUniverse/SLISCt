# Makefile

exe = slisc.x
eigenPath = ../EigenTest/Eigen/

# this project requires "MatFile_linux" github repository

source = \
main.cpp \
disp.cpp \
fft.cpp \
interp1.cpp \
eigen_basics.cpp \
eigen_fft.cpp \
eigen_linsolve.cpp

objects = $(source:.cpp=.o)

compiler = g++

flags =  \
-I $(eigenPath) \
-std=c++11 \
-O3
# -g
# -O3

$(exe):$(objects)
	$(compiler) -o $(exe) $(flags) $(objects) $(xobjects) $(libs)

$(objects):$(source)
	$(compiler) -c $(flags) $(source) $(sourcepath)$(xsource)

clean:
	rm -f *.o *.x *.gch
