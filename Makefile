# Makefile

exe = main.x

eigenPath = ../EigenTest/Eigen/

Sdir = SLISC/gsl/source/

# this project requires "MatFile_linux" github repository

source = main.cpp #SLISC/print.cpp

compiler = icpc
# (NOTE:  The  icpc command uses the same compiler options as the icc com-
# mand. Invoking the compiler using icpc compiles .c and .i files as C++.
# Invoking  the  compiler  using icc compiles .c and .i files as C. Using
# icpc always links in  C++  libraries.  Using  icc  only  links  in  C++
# libraries if C++ source is provided on the command line.)

flags = -Wall -I $(eigenPath) -I SLISC/gsl/include -I SLISC/gsl/source -std=c++17 -g -mkl -fp-model precise -fp-model except
# -O3

$(exe): main.o
	$(compiler) -o $(exe) $(flags) *.o

main.o: main.cpp coulomb.o SLISC/*.h test/*.h
	$(compiler) -c $(flags) $(source)

clean:
	rm -f *.o *.x *.gch

stream.o: $(Sdir)err/stream.c
	$(compiler) $(flags) -c $<
error.o: $(Sdir)err/error.c stream.o
	$(compiler) $(flags) -c $<
log.o: $(Sdir)specfunc/log.c  error.o
	$(compiler) $(flags) -c $<
pow_int.o: pow_int.c error.o
	$(compiler) $(flags) -c $<
trig.o: $(Sdir)specfunc/trig.c log.o error.o
	$(compiler) $(flags) -c $<
math.o: $(Sdir)complex/math.c
	$(compiler) $(flags) -c $<
fdiv.o: $(Sdir)sys/fdiv.c
	$(compiler) $(flags) -c $<
infnan.o: $(Sdir)sys/infnan.c fdiv.o
	$(compiler) $(flags) -c $<
coerce.o: $(Sdir)sys/coerce.c
	$(compiler) $(flags) -c $<
elementary.o: $(Sdir)specfunc/elementary.c coerce.o error.o
	$(compiler) $(flags) -c $<
exp.o: $(Sdir)specfunc/exp.c error.o gamma.o
	$(compiler) $(flags) -c $<
zeta.o: $(Sdir)specfunc/zeta.c exp.o error.o gamma.o elementary.o
	$(compiler) $(flags) -c $<
psi.o: $(Sdir)specfunc/psi.c math.o exp.o error.o infnan.o gamma.o zeta.o
	$(compiler) $(flags) -c $<
gamma.o: $(Sdir)specfunc/gamma.c  log.o exp.o error.o psi.o trig.o
	$(compiler) $(flags) -c $<
poch.o: poch.c log.o exp.o error.o gamma.o pow_int.o psi.o
	$(compiler) $(flags) -c $<
bessel_temme.o: bessel_temme.c
	$(compiler) $(flags) -c $<
bessel_amp_phase.o: bessel_amp_phase.c
	$(compiler) $(flags) -c $<
bessel.o: bessel.c exp.o error.o poch.o bessel_amp_phase.o gamma.o elementary.o bessel_temme.o
	$(compiler) $(flags) -c $<
bessel_J1.o: bessel_J1.c  error.o bessel.o
	$(compiler) $(flags) -c $<
airy_der.o: airy_der.c exp.o error.o
	$(compiler) $(flags) -c $<
airy.o: $(Sdir)specfunc/airy.c error.o trig.o
	$(compiler) $(flags) -c $<
bessel_olver.o: bessel_olver.c  airy_der.o airy.o error.o
	$(compiler) $(flags) -c $<
bessel_J0.o: bessel_J0.c  error.o bessel.o
	$(compiler) $(flags) -c $<
bessel_Jn.o: bessel_Jn.c bessel_J1.o error.o bessel.o bessel_olver.o bessel_J0.o
	$(compiler) $(flags) -c $<
legendre_P.o: $(Sdir)specfunc/legendre_P.c
	$(compiler) $(flags) -c $<
legendre_poly.o: legendre_poly.c log.o error.o poch.o bessel_Jn.o bessel_J0.o
	$(compiler) $(flags) -c $<
laguerre.o: laguerre.c error.o exp.o gamma.o
	$(compiler) $(flags) -c $<
coupling.o: $(Sdir)specfunc/coupling.c exp.o gamma.o error.o infnan.o
	$(compiler) $(flags) -c $<
coulomb_bound.o: coulomb_bound.c laguerre.o pow_int.o gamma.o
	$(compiler) $(flags) -c $<
coulomb.o: $(Sdir)specfunc/coulomb.c exp.o airy.o psi.o gamma.o error.o infnan.o
	$(compiler) $(flags) -c $<
