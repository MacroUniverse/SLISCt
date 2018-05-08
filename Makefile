#
# Makefile
#
# Author: Uwe Thumm
# Date:   6/5/99
#




FOBJS =   \
main.o

COMPILER = g++

FFLAGS = -g
#FFLAGS = -O3

.cpp.o:
	$(COMPILER) -c $(FFLAGS) *.cpp *.h


nr3: $(FOBJS) $(LIB)
	$(COMPILER) $(FFLAGS) -o nr3.x $(FOBJS)

clean:
	rm -f *.o *.x *.gch
