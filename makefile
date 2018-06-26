##############################################################################
#  Makefile for libmetreader.a
#
#    User-specified flags are in this top block
#
###############################################################################

#      This file is a component of the volcanic ash transport and dispersion model Ash3d,
#      written at the U.S. Geological Survey by Hans F. Schwaiger (hschwaiger@usgs.gov),
#      Larry G. Mastin (lgmastin@usgs.gov), and Roger P. Denlinger (roger@usgs.gov).

#      The model and its source code are products of the U.S. Federal Government and therefore
#      bear no copyright.  They may be copied, redistributed and freely incorporated 
#      into derivative products.  However as a matter of scientific courtesy we ask that
#      you credit the authors and cite published documentation of this model (below) when
#      publishing or distributing derivative products.

#      Schwaiger, H.F., Denlinger, R.P., and Mastin, L.G., 2012, Ash3d, a finite-
#         volume, conservative numerical model for ash transport and tephra deposition,
#         Journal of Geophysical Research, 117, B04204, doi:10.1029/2011JB008968. 

#      We make no guarantees, expressed or implied, as to the usefulness of the software
#      and its documentation for any purpose.  We assume no responsibility to provide
#      technical support to users of this software.

#  SYSTEM specifies which compiler to use
#    Current available options are:
#      gfortran , ifort
#    This variable cannot be left blank
#      
SYSTEM = gfortran
#
#  RUN specifies which collection of compilation flags that should be run
#    Current available options are:
#      DEBUG : includes debugging info and issues warnings
#      PROF  : includes profiling flags with some optimization
#      OPT   : includes optimizations flags for fastest runtime
#    This variable cannot be left blank
RUN = OPT
#
INSTALLDIR=/opt/USGS
#INSTALLDIR=~/gcc
#
# DATA FORMATS
#  For each data format you want to include in the library, set the corresponding
#  variable below to 'T'.  Set to 'F' any you do not want compiled or any unavailable
USENETCDF = T
USEGRIB2 = T
USEHDF = F   # Note: the hdf reader is not yet functional

# MEMORY
# If you need pointer arrays instead of allocatable arrays, set this to 'T'
USEPOINTERS = F

###############################################################################
#####  END OF USER SPECIFIED FLAGS  ###########################################
###############################################################################

FPPFLAGS = 
ifeq ($(USENETCDF), T)
 ncFPPFLAG = -DUSENETCDF
 ncOBJS = MetReader_NetCDF.o
 nclib = -lnetcdf -lnetcdff
else
 ncFPPFLAG =
 ncOBJS =
 nclib =
endif
ifeq ($(USEGRIB2), T)
 grb2FPPFLAG = -DUSEGRIB2
 grb2OBJS = MetReader_GRIB.o MetReader_GRIB_index.o
 grblib = -lgrib_api_f90 -lgrib_api
else
 grb2FPPFLAG =
 grb2OBJS =
 grblib =
endif
ifeq ($(USEHDF), T)
 hdfFPPFLAG = -DUSEHDF
else
 hdfFPPFLAG =
endif

ifeq ($(USEPOINTERS), T)
 memFPPFLAG = -DUSEPOINTERS
else
 memFPPFLAG =
endif

FPPFLAGS = -x f95-cpp-input $(ncFPPFLAG) $(grb2FPPFLAG) $(hdfFPPFLAG) $(memFPPFLAG)

# location of HoursSince and projection
USGSLIBDIR = -L$(INSTALLDIR)/lib
USGSINC = -I$(INSTALLDIR)/include
USGSLIB = $(USGSLIBDIR) $(USGSINC) -lhourssince -lprojection

EXEC = \
 gen_GRIB2_index   \
 tools/ncMetSonde  \
 tools/ncMetTraj_F \
 tools/ncMetTraj_B \
 tools/ncMet_check

###############################################################################
###############################################################################

###############################################################################
##########  GNU Fortran Compiler  #############################################
ifeq ($(SYSTEM), gfortran)
    FCHOME=/usr
    FC = /usr/bin/gfortran

    COMPINC = -I$(FCHOME)/include -I$(FCHOME)/lib64/gfortran/modules -I$(INSTALLDIR)/include
    COMPLIBS = -L$(FCHOME)/lib64 -L${INSTALLDIR}/lib 

    LIBS = $(COMPLIBS) $(COMPINC)
    # -lefence 
# Debugging flags
ifeq ($(RUN), DEBUG)
    FFLAGS =  -O0 -g3 -Wall -fbounds-check -pedantic -fimplicit-none -Wunderflow -Wuninitialized -ffpe-trap=invalid,zero,overflow -fdefault-real-8 
endif
# Profiling flags
ifeq ($(RUN), PROF)
    FFLAGS = -g -pg -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
# Production run flags
ifeq ($(RUN), OPT)
    FFLAGS = -O3 -w -fno-math-errno -funsafe-math-optimizations -fno-trapping-math -fno-signaling-nans -fcx-limited-range -fno-rounding-math -fdefault-real-8
endif
    EXFLAGS =
endif
###############################################################################

all: libMetReader.a tools

libMetReader.a: MetReader.F90 MetReader.o $(ncOBJS) $(grb2OBJS) MetReader_Grids.o MetReader_ASCII.o makefile
	ar rcs libMetReader.a MetReader.o $(ncOBJS) $(grb2OBJS) MetReader_Grids.o MetReader_ASCII.o

MetReader.o: MetReader.F90 makefile
	$(FC) $(FPPFLAGS) $(EXFLAGS) -c MetReader.F90
MetReader_Grids.o: MetReader_Grids.f90 MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c MetReader_Grids.f90
MetReader_ASCII.o: MetReader_ASCII.f90 MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -c MetReader_ASCII.f90

ifeq ($(USENETCDF), T)
MetReader_NetCDF.o: MetReader_NetCDF.f90 MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) -c MetReader_NetCDF.f90
endif
ifeq ($(USEGRIB2), T)
MetReader_GRIB_index.o: MetReader_GRIB_index.f90 makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c MetReader_GRIB_index.f90
MetReader_GRIB.o: MetReader_GRIB.f90 MetReader_GRIB_index.o MetReader.o makefile
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c MetReader_GRIB.f90
gen_GRIB2_index: gen_GRIB2_index.f90 MetReader_GRIB_index.o makefile libMetReader.a
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(grblib) -c gen_GRIB2_index.f90
	$(FC) $(FFLAGS) $(EXFLAGS) MetReader_GRIB_index.o gen_GRIB2_index.o $(LIBS) $(grblib) -o gen_GRIB2_index
endif
#ifeq ($(USEHDF), T)
#MetReader_HDF.o: MetReader_HDF.f90 MetReader.o makefile
#	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) -lhdf -Iinclude/ -c MetReader_HDF.f90
#endif


ifeq ($(USEGRIB2), T)
  GRIBTOOL = gen_GRIB2_index
else
  GRIBTOOL =
endif

tools: ncMetSonde ncMetTraj_F ncMetTraj_B ncMet_check $(GRIBTOOL)

ncMetSonde: tools/ncMetSonde.f90 makefile libMetReader.a
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/ncMetSonde.f90
	$(FC) $(FFLAGS) $(EXFLAGS) ncMetSonde.o $(LIBS) $(nclib) $(grblib) -L./ -lMetReader $(USGSLIB) -o tools/ncMetSonde
ncMetTraj_F: tools/ncMetTraj.F90 makefile libMetReader.a
	$(FC) -x f95-cpp-input -DFORWARD  $(FFLAGS) $(EXFLAGS) tools/ncMetTraj.F90 -o tools/ncMetTraj_F $(LIBS) $(nclib) $(grblib) -L./ -lMetReader $(USGSLIB)
ncMetTraj_B: tools/ncMetTraj.F90 makefile libMetReader.a
	$(FC) -x f95-cpp-input -DBACKWARD $(FFLAGS) $(EXFLAGS) tools/ncMetTraj.F90 -o tools/ncMetTraj_B $(LIBS) $(nclib) $(grblib) -L./ -lMetReader $(USGSLIB)
ncMet_check: tools/ncMet_check.f90 makefile libMetReader.a
	$(FC) $(FFLAGS) $(EXFLAGS) $(LIBS) $(nclib) $(grblib) -c tools/ncMet_check.f90
	$(FC) $(FFLAGS) $(EXFLAGS) ncMet_check.o $(LIBS) $(nclib) $(grblib) -L./ -lMetReader $(USGSLIB) -o tools/ncMet_check

clean:
	rm -f *.o
	rm -f *.mod
	rm -f lib*.a
	rm -f $(EXEC)

install:
	install -d $(INSTALLDIR)/lib/
	install -d $(INSTALLDIR)/include/
	install -d $(INSTALLDIR)/bin/
	install -d $(INSTALLDIR)/bin/autorun_scripts
	install -m 644 libMetReader.a $(INSTALLDIR)/lib/
	install -m 644 *.mod $(INSTALLDIR)/include/
	install -m 755 $(EXEC) $(INSTALLDIR)/bin/
	install -m 755 autorun_scripts/*.sh $(INSTALLDIR)/bin/autorun_scripts/
