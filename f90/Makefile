# Compiler:
F90C  = mpif90
CC  = mpicc


HEALPIX=$(HEALPIX3)



# Link libraries and options
FFLAGS  = -O3 -openmp -openmp_report0 -fPIC -I$(HEALPIX)/include -I$(HEALPIX)/src/f90/mod 
FFLAGS_PAR = -03 -heap-arrays 16 -I$(HEALPIX)/include  -I$(CFITSIO)/include 
CFLAGS	=-O2 -fomit-frame-pointer -fforce-addr -fforce-mem
LDFLAGS	= -L$(HEALPIX)/lib -L$(CFITSIO)/lib -lhealpix -lgif -lcfitsio -L/opt/app-sync/intel/composer_xe_2013_sp1.2.144/compiler/lib/intel64/  -limf


#

#LDFLAGS = -lmpi -lfmpi 
INTLIB = -L/mn/stornext/u2/yabebalf/lib -lintlib 


# Set up appropriate suffix list
.SUFFIXES:
.SUFFIXES: .o .f .for .f90 .c

# Compile-time options
# F90 options

BIN=../bin


# Obj list:
OBJ	=  modif_pixtools.o
OBJW =  get_wmap_dipdir.o


OBJPIX = modif_pixtools.o
#+
# Make some default rules
#-

# Double suffix rules:
.f90.o .f.o .for.o :
	$(F90C) $(FFLAGS) -c $<
.c.o	:
	$(CC) $(CFLAGS) -c $<

# Single suffix rules:
.f90 .f .for :
	$(F90C) $(FFLAGS) -o $@ $<

.o :
	$(F90C) $(LDFLAGS) -o $@ $<



# Build all programs
all: get_dipdir get_wmap_dipdir get_cross_dipdir get_cross_lmin_dipdir get_cross_lmax_dipdir get_listpix_disk


query_discvec_test: ${OBJ}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${OBJ} ${LDFLAGS}


npixels_per_disk: ${OBJ}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${OBJ} ${LDFLAGS}

get_listpix_ring_disk: ${OBJ}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${OBJ} ${LDFLAGS}

get_listpix_nest_disk: ${OBJ}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${OBJ} ${LDFLAGS}

#---------

get_dipdir: 
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}

get_wmap_dipdir: ${OBJW}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}

get_cross_dipdir: ${OBJW}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}

get_variance_dipdir: ${OBJW}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}

get_cross_lmin_dipdir: ${OBJW}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}

get_cross_lmax_dipdir: ${OBJW}
		$(F90C) $(FFLAGS) -o $(HEALPIX)/bin/$@ $@.f90 ${LDFLAGS}
# Remove all files created from building 
clean	:
	-$(RM) hess
	-$(RM) *.mod
	-$(RM) *.o
	-$(RM) *.lst

tidy:	clean

