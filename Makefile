obj =  wmpi.o module_4sites.o f90code_spinless_allo.o 

# compiler
# Here mpif90 should be compiled with gfortran
f90  = mpif90 -cpp -DMPI -DINTELMKL -g

flag = -O3                              

# blas and lapack libraries
libs = -L ${OPENBLAS_ROOT}/lib -lopenblas\

 
main :  $(obj)
	$(f90) $(obj) -o exactdiag.x $(libs) 
	cp -f exactdiag.x ../bin

.SUFFIXES: .o .f90

.f90.o :
	$(f90) -c $(flag) $(includes) $*.f90

clean :
	rm -f *.o *.mod *~ *.x
