# Makefile

exe = nr3.x


# this project requires "MatFile_linux" github repository

source = \
main.cpp \
nr3plus.cpp \
interp_1d.cpp \
interp_2d.cpp \
tridag.cpp

objects = $(source:.cpp=.o)

compiler = g++

flags =  \
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

