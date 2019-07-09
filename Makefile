# Makefile

exe = main.x

eigenPath = ../EigenTest/Eigen/

Sdir = SLISC/gsl/source/

# this project requires "MatFile_linux" github repository

compiler = icpc
# (NOTE:  The  icpc command uses the same compiler options as the icc com-
# mand. Invoking the compiler using icpc compiles .c and .i files as C++.
# Invoking  the  compiler  using icc compiles .c and .i files as C. Using
# icpc always links in  C++  libraries.  Using  icc  only  links  in  C++
# libraries if C++ source is provided on the command line.)

flags = -Wall -I $(eigenPath) -I SLISC/gsl/include -I SLISC/gsl/source -std=c++17 -g -mkl -fp-model precise -fp-model except -qopenmp
# -O3 # highest optimization
# -qopenmp # run OpenMP in parallel mode
# -qopenmp-stubs # run OpenMP in serial mode
# -fp-model # floating point model

# goal
$(exe): main.o
	$(compiler) -o $(exe) $(flags) *.o

clean:
	rm -f *.o *.x *.gch

MAKEFLAGS = -r # no default implicit rule
.SUFFIXES: .c .cpp .o .x

# implicit rules (paths must be specific)
# there can be more dependencies than specified here

%.o: %.cpp
	$(compiler) $(flags) -c $<

%.o: $(Sdir)err/%.c
	$(compiler) $(flags) -c $<

%.o: $(Sdir)specfunc/%.c
	$(compiler) $(flags) -c $<

%.o: $(Sdir)complex/%.c
	$(compiler) $(flags) -c $<

%.o: $(Sdir)sys/%.c
	$(compiler) $(flags) -c $<

# main file dependency (must include any used GSL function's object file)
main.o: main.cpp SLISC/*.h test/*.h coulomb.o coupling.o coulomb_bound.o legendre_poly.o legendre_P.o

# GSL source dependency (use Understand cluster graph to add)
# the object files on the right are not needed to compile, but needed to link
# this way, only needed object files are compiled, no matter how long the list below is

stream.o: $(Sdir)err/stream.c
error.o: $(Sdir)err/error.c stream.o
log.o: $(Sdir)specfunc/log.c  error.o
pow_int.o: $(Sdir)sys/pow_int.c error.o
trig.o: $(Sdir)specfunc/trig.c log.o error.o
math.o: $(Sdir)complex/math.c
fdiv.o: $(Sdir)sys/fdiv.c
infnan.o: $(Sdir)sys/infnan.c fdiv.o
coerce.o: $(Sdir)sys/coerce.c
elementary.o: $(Sdir)specfunc/elementary.c coerce.o error.o
exp.o: $(Sdir)specfunc/exp.c error.o gamma.o
zeta.o: $(Sdir)specfunc/zeta.c exp.o error.o gamma.o elementary.o
psi.o: $(Sdir)specfunc/psi.c math.o exp.o error.o infnan.o gamma.o zeta.o
gamma.o: $(Sdir)specfunc/gamma.c  log.o exp.o error.o psi.o trig.o
poch.o: $(Sdir)specfunc/poch.c log.o exp.o error.o gamma.o pow_int.o psi.o
bessel_temme.o: $(Sdir)specfunc/bessel_temme.c
bessel_amp_phase.o: $(Sdir)specfunc/bessel_amp_phase.c
bessel.o: $(Sdir)specfunc/bessel.c exp.o error.o poch.o bessel_amp_phase.o gamma.o elementary.o bessel_temme.o
bessel_J1.o: $(Sdir)specfunc/bessel_J1.c error.o bessel.o
airy_der.o: $(Sdir)specfunc/airy_der.c exp.o error.o
airy.o: $(Sdir)specfunc/airy.c error.o trig.o
bessel_olver.o: $(Sdir)specfunc/bessel_olver.c  airy_der.o airy.o error.o
bessel_J0.o: $(Sdir)specfunc/bessel_J0.c error.o bessel.o
bessel_Jn.o: $(Sdir)specfunc/bessel_Jn.c bessel_J1.o error.o bessel.o bessel_olver.o bessel_J0.o
legendre_P.o: $(Sdir)specfunc/legendre_P.c
legendre_poly.o: $(Sdir)specfunc/legendre_poly.c log.o error.o poch.o bessel_Jn.o bessel_J0.o
laguerre.o: $(Sdir)specfunc/laguerre.c error.o exp.o gamma.o
coupling.o: $(Sdir)specfunc/coupling.c exp.o gamma.o error.o infnan.o
coulomb_bound.o: $(Sdir)specfunc/coulomb_bound.c laguerre.o pow_int.o gamma.o
coulomb.o: $(Sdir)specfunc/coulomb.c exp.o airy.o psi.o gamma.o error.o infnan.o
