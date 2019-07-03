# Makefile

exe = main.x

eigenPath = ../EigenTest/Eigen/

Sdir = gsl/source/

# this project requires "MatFile_linux" github repository

source = main.cpp #SLISC/print.cpp
objects = main.o #print.o
gslobjects = math.o error.o stream.o coupling.o elementary.o exp.o gamma.o log.o psi.o trig.o zeta.o legendre_P.o coerce.o fdiv.o infnan.o

compiler = icpc
# (NOTE:  The  icpc command uses the same compiler options as the icc com-
# mand. Invoking the compiler using icpc compiles .c and .i files as C++.
# Invoking  the  compiler  using icc compiles .c and .i files as C. Using
# icpc always links in  C++  libraries.  Using  icc  only  links  in  C++
# libraries if C++ source is provided on the command line.)

flags = -Wall -I $(eigenPath) -I gsl/include -I gsl/source -std=c++17 -g -mkl -fp-model precise -fp-model except
# -O3

$(exe):$(objects) $(gslobjects)
	$(compiler) -o $(exe) $(flags) $(objects) $(gslobjects)

$(objects):$(source)
	$(compiler) -c $(flags) $(source)

clean:
	rm -f *.o *.x *.gch

math.o: $(Sdir)complex/math.c
	$(compiler) $(flags) -c $<

error.o: $(Sdir)err/error.c
	$(compiler) $(flags) -c $<
stream.o: $(Sdir)err/stream.c
	$(compiler) $(flags) -c $<

coupling.o: $(Sdir)specfunc/coupling.c
	$(compiler) $(flags) -c $<
elementary.o: $(Sdir)specfunc/elementary.c
	$(compiler) $(flags) -c $<
exp.o: $(Sdir)specfunc/exp.c
	$(compiler) $(flags) -c $<
gamma.o: $(Sdir)specfunc/gamma.c
	$(compiler) $(flags) -c $<
log.o: $(Sdir)specfunc/log.c
	$(compiler) $(flags) -c $<
psi.o: $(Sdir)specfunc/psi.c
	$(compiler) $(flags) -c $<
trig.o: $(Sdir)specfunc/trig.c
	$(compiler) $(flags) -c $<
zeta.o: $(Sdir)specfunc/zeta.c
	$(compiler) $(flags) -c $<
legendre_P.o: $(Sdir)specfunc/legendre_P.c
	$(compiler) $(flags) -c $<

coerce.o: $(Sdir)sys/coerce.c
	$(compiler) $(flags) -c $<
fdiv.o: $(Sdir)sys/fdiv.c
	$(compiler) $(flags) -c $<
infnan.o: $(Sdir)sys/infnan.c
	$(compiler) $(flags) -c $<
