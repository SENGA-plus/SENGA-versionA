#
# make file for SENGA2
#
# RSC 08-JAN-2005 LAM/MPI implementation
# RSC 12-JUN-2008 SUSE/MPICH implementation
#

# Fortran compiler
#FC = mpif90# -fallow-argument-mismatch -g -fcheck=all -cpp
FC = h5pfc -g -fcheck=all -cpp#-DHDF5 

# Quick non-optimized compile
#FFLAGS = -g -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Non-optimized compile with profiling
#FFLAGS = -pg -fno-silent -ftrapping-math -ftrapv -Wimplicit -Wall

# Optimized compile
FFLAGS = -O3 -Wall -mcmodel=medium

# Optimized compile with profiling
#FFLAGS = -pg -O3 -fno-silent -Wimplicit -Wall

# Source files
SRC =   senga2.f rhscal.f rhsvel.f lincom.f fincom.f \
	temper.f tempin.f chrate.f chemin.f diffin.f radcal.f radcin.f \
	hdf5io.f indata.f contin.f output.f finish.f \
	dfinit.f dfbydx.f dfbydy.f dfbydz.f d2fdx2.f d2fdy2.f d2fdz2.f \
	d2fdxy.f d2fdxz.f d2fdyz.f dfmstr.f \
	flinit.f fltrem.f filtrx.f filtry.f filtrz.f \
	dtinit.f adaptt.f \
	turbin.f turbft.f fftsym.f fftf3d.f \
	espect.f espini.f espovk.f espksq.f integf.f \
	flamin.f erfunc.f \
	bcinit.f bounds.f boundt.f bountt.f \
	bcutxl.f bcttxl.f bcdtxl.f bcytxl.f bctixl.f buftxl.f fftixl.f \
	bcutxr.f bcttxr.f bcdtxr.f bcytxr.f \
	bcutyl.f bcttyl.f bcdtyl.f bcytyl.f \
	bcutyr.f bcttyr.f bcdtyr.f bcytyr.f \
	bcutzl.f bcttzl.f bcdtzl.f bcytzl.f \
	bcutzr.f bcttzr.f bcdtzr.f bcytzr.f \
	pardom.f parfer.f \
	cxhalo.f cyhalo.f czhalo.f cxhals.f cyhals.f czhals.f \
	exhalo.f inhalo.f exhals.f inhals.f \
	zeroxl.f zeroxr.f zeroyl.f zeroyr.f zerozl.f zerozr.f \
	ardump.f ardums.f chkarr.f chkary.f

# Include files
INC = com_senga2.h com_espect.h

# Libraries: need separate compilation in their own subdirectories
# FFTs; random number generators
# MPI libraries
LIBS = fftlib/fftlib.a rndlib/rndlib.a parlib/parlib.a -lm

# Object files are the same as the source files but with .o rather than .f
OBJS=$(SRC:.f=.o)

# Default compilation of .f files to .o
.f.o: $(FC) $(INC) -c $<

# Full make of SENGA2
all: senga2

senga2: $(OBJS) 
	$(FC) $(OBJS) -o senga2 $(LIBS)

# Add -pg for profiling
#senga2: $(OBJS) 
#	$(FC) $(OBJS) -pg -o senga2 $(LIBS)

# Command to clean up after
clean: 
	rm -f *.o senga2 core *.mod
