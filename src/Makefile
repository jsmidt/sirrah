#
#  Compiler parameters
#
FC	= gfortran
FFLAGS	= -O2

#
# fortranlib [1] by astrofrog is needed to compile this code. Put the
# path to this code in FLIBDIR. Do not edit FLIB and FOBS.
# [1] https://github.com/astrofrog/fortranlib
#
FLIBDIR = ../../fortranlib
FLIB    = $(FLIBDIR)/src/base_types.o $(FLIBDIR)/src/lib_constants.o $(FLIBDIR)/src/lib_array.o
FOBS = base_types.o lib_constants.o lib_array.o

# The files needed to compile unique to SIRRAH.
OBJS	= io.o hydro_1d.o driver.o 	

# The compile tree
default: sirrah 
all: sirrah

sirrah: $(FLIB) $(OBJS) 
	$(FC) $(FFLAGS) -o ../sirrah $(OBJS) $(FOBS)
	
# Compile instructions
.SUFFIXES: .F .f90
.F.o: 
	$(FC) $(FFLAGS) -c $<
.f90.o: 
	$(FC) $(FFLAGS) -c $<

# Clean
clean:
	-rm -f *.o *~ *.mod ../sirrah
