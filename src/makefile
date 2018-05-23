###################################################################
###################################################################
#
# Basic makefile to compile the En-GARD circulation downscaling code.
#
# To compile this code on your machine, edit the NCDF_PATH variable and the fortran compiler
#
# To build the model simply run make.
#
# To specify the location of the NETCDF or LAPACK libraries, or the fortran compiler
#   make NETCDF=/usr/local LAPACK=/usr/local FC=gfortran
#
# Parallel compiles can speed up the compile time considerably, consider
# 	make -j4
#
# The primary user facing rules are:
# 	all		: builds the model, this is the default
# 	clean	: remove any build objects (not executables)
#  	allclean: make clean and remove executables
# 	install	: make the executable and copy it to ~/bin (or other install dir)
# 	test	: build various test executables
#
###################################################################
###################################################################

ifeq ($(FC),)
	FC=gfortran-6
endif

ifeq ($(COMPILER),)
	ifeq ($(FC), gfortran)
		COMPILER=gnu
	endif
	ifeq ($(FC), gfortran-6)
		COMPILER=gnu
	endif
	ifeq ($(FC), /usr/bin/gfortran-6)
		COMPILER=gnu
	endif
	ifeq ($(FC), ifort)
		COMPILER=intel
	endif
	ifeq ($(FC), pgf90)
		COMPILER=pgi
	endif
	ifeq ($(COMPILER), cray)
		COMPILER=cray
	endif
endif

ifeq ($(NETCDF),)
	NETCDF = $(shell nc-config --prefix)
endif
ifeq ($(NETCDF),)
	NETCDF = /usr/local
endif
NCDF_PATH = $(NETCDF)

ifeq ($(LAPACK),)
	LAPACK=/usr/local
endif
LAPACK_PATH = $(LAPACK)

INSTALL_DIR = ~/bin
RM=/bin/rm
CP=/bin/cp
BUILD=build/

###################################################################
###################################################################
#
#	Try to identify the computer we are running on and set up library paths
#
###################################################################
###################################################################
NODENAME := $(shell uname -n)

ifeq ($(patsubst cheyenne%,cheyenne,$(NODENAME)), cheyenne)
	NODENAME=cheyenne
endif

ifeq ($(NODENAME), cheyenne)
	# on cheyenne:
    ifeq ($(LMOD_FAMILY_COMPILER),gnu)
        FC=gfortran
		COMPILER=gnu
        NCDF_PATH=/glade/u/apps/ch/opt/netcdf/4.4.1.1/gnu/6.2.0
        LAPACK_PATH=/glade/u/apps/ch/opt/netlib/3.7.0/gnu/4.8.5
    endif
	ifeq ($(LMOD_FAMILY_COMPILER),intel)
	    FC=ifort
		COMPILER=intel
	    NCDF_PATH=/glade/u/apps/ch/opt/netcdf/4.4.1.1/intel/17.0.1
		LAPACK_PATH=/glade/u/apps/opt/intel/2017u1/compilers_and_libraries/linux/mkl/lib ${NCAR_LIBS_MKL} \#
		# this is an ugly hack to get the NCAR_LIBS_MKL variable in, need to rethink a little
    endif
endif


# get GIT version info
GIT_VERSION := $(shell git describe --long --dirty --all --always | sed -e's/heads\///')
UNAME_S := $(shell uname -s)
HOSTNAME := $(shell uname -n)

###################################################################
###################################################################
#
#	Compiler specific compilation flags
#
###################################################################
###################################################################
ifeq ($(COMPILER), gnu)
	COMP=-fopenmp -lgomp -O3 -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops \
		-fno-protect-parens -ffast-math #-flto #
	LINK=-fopenmp -lgomp
	PREPROC=-cpp
	MODOUTPUT=-J $(BUILD)
endif
ifeq ($(COMPILER), gnu)
	COMP=-fopenmp -lgomp -O3 -ffree-line-length-none -ftree-vectorize -fimplicit-none -funroll-loops \
		-fno-protect-parens -ffast-math #-flto #
	LINK=-fopenmp -lgomp
	PREPROC=-cpp
	MODOUTPUT=-J $(BUILD)
endif
# Intel fortran
ifeq ($(COMPILER), intel)
	COMP= -u -qopenmp -O3 -no-prec-div -xHost -ftz
	LINK= -liomp5
	PREPROC=-fpp
	MODOUTPUT=-module $(BUILD)
endif
# PGI fortran
ifeq ($(COMPILER), pgi)
	COMP=-O2 -mp -Mdclchk #-fast -O3 -mp -c -Mdclchk
	LINK=-mp
	PREPROC=-Mpreprocess
	MODOUTPUT=-module $(BUILD)
endif

# Cray fortran
ifeq ($(COMPILER), cray)
	COMP= -h omp vector2 -O2 -eI
	LINK= -fopenmp
	PREPROC= -eZ
	MODOUTPUT= -J $(BUILD) -em
endif

###################################################################
###################################################################
#
# Various compiling options.  Set the MODE variable with "make MODE=debugslow" etc.
#
###################################################################
###################################################################
ifeq ($(MODE), debugslow)
	ifeq ($(COMPILER), intel)
		COMP= -debug -debug-parameters all -traceback -ftrapuv -g -fpe0 -u -check all -check noarg_temp_created -CB
		LINK=
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -g -Wall -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK=
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -g -Wall -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK=
	endif
	ifeq ($(COMPILER), pgi)
		COMP= -g -Mbounds -Mlist -Minfo  -Mdclchk
		LINK=
	endif
	ifeq ($(COMPILER), cray)
		COMP=-h noomp -g -m 0 -R abcsp
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debug)
	ifeq ($(COMPILER), intel)
		COMP= -debug -O1 -u -check all -check noarg_temp_created -traceback -fpe0 -fast-transcendentals -xhost
		LINK=
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -O1 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK=
	endif
	ifeq ($(COMPILER), pgi)
		COMP= -gopt -O1 -Mbounds -Mlist -Minfo  -Mdclchk
		LINK=
	endif
	ifeq ($(COMPILER), cray)
		COMP=-O1 -h noomp -g
		LINK=-h noomp
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugompslow)
	ifeq ($(COMPILER), intel)
		COMP= -qopenmp -debug -u	-fpe0 -traceback -check all -check noarg_temp_created -fp-stack-check
		LINK= -liomp5
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -fopenmp -lgomp -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK= -fopenmp -lgomp
	endif
	ifeq ($(COMPILER), pgi)
		COMP= -g -Mbounds -Mlist -Minfo -mp -Mdclchk
		LINK= -mp
	endif
	ifeq ($(COMPILER), cray)
		COMP= -g -m 0 -R abcsp
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif
ifeq ($(MODE), debugomp)
	ifeq ($(COMPILER), intel)
		COMP= -qopenmp -debug -O3 -u -traceback -fpe0 -ftz -fast-transcendentals  -xHost
		LINK= -liomp5
	endif
	ifeq ($(COMPILER), gnu)
		COMP= -fopenmp -lgomp -O1 -g -fbounds-check -fbacktrace -finit-real=nan -ffree-line-length-none
		LINK= -fopenmp -lgomp
	endif
	ifeq ($(COMPILER), pgi)
		COMP= -g -O1 -Mbounds -Mlist -Minfo -mp -Mdclchk
		LINK= -mp
	endif
	ifeq ($(COMPILER), cray)
		COMP=-O1 -g
		LINK=
		PREPROC=-eZ
		MODOUTPUT=-e m -J $(BUILD)
	endif
endif

PROF=
ifeq ($(MODE), profile)
	ifeq ($(COMPILER), intel)
		PROF=-pg -debug inline-debug-info -shared-intel
		COMP= -u -qopenmp -O3 -no-prec-div -xHost -ftz #because -fast includes -static
		# not available in ifort <13 -align array64byte
	endif
	ifeq ($(COMPILER), gnu)
		PROF=-g
	endif
endif
ifeq ($(MODE), fast) # WARNING -ipo (included in -fast) takes forever for very little gain, and this may be unstable
	ifeq ($(COMPILER), intel)
		COMP= -u -qopenmp -fast -ftz #-fast-transcendentals # not available in ifort <13: -align array64byte
	endif
	ifeq ($(COMPILER), gnu)
		COMP=-O3 -finit-real=nan -ffree-line-length-none
		LINK=
	endif

endif

$(info $$NODENAME    = ${NODENAME})
$(info $$FC          = ${FC})
$(info $$NCDF_PATH   = ${NCDF_PATH})
$(info $$LAPACK_PATH = ${LAPACK_PATH})
$(info $$GIT_VERSION = ${GIT_VERSION})
$(info $$COMP        = ${COMP})
$(info $$LINK        = ${LINK})
$(info $$PREPROC     = ${PREPROC})
$(info $$MODOUTPUT   = ${MODOUTPUT})
$(info $$MODE        = ${MODE})

###################################################################
###################################################################
#
# Should not need to edit anything below this line
#
###################################################################
###################################################################
# copy required libraries into a directory accessible on compute nodes and set LD_RUN_PATH e.g.
# export LD_RUN_PATH=$LD_RUN_PATH:/path/to/netcdf/libraries/lib:/path/to/fftw/libraries/lib

# flags forall
FCFLAGS = $(COMP) $(PROF) -I$(NCDF_PATH)/include -I$(BUILD) $(INC_NCAR)

# libraries needed for linking
LDFLAGS = $(LINK) $(PROF) -L$(NCDF_PATH)/lib -lnetcdf -lnetcdff $(LIB_NCAR) -L$(LAPACK_PATH)/lib -llapack #-lblas

OBJS =  \
	$(BUILD)time.o            	\
	$(BUILD)string.o          	\
	$(BUILD)model_constants.o 	\
	$(BUILD)data_structures.o 	\
	$(BUILD)configuration.o   	\
	$(BUILD)init.o			  	\
	$(BUILD)io_routines.o     	\
	$(BUILD)io_routines_h.o    	\
	$(BUILD)time_io.o         	\
	$(BUILD)time_utils.o       	\
	$(BUILD)atmosphere_io.o	  	\
	$(BUILD)atmosphere_io_h.o	\
	$(BUILD)gcm_io.o          	\
	$(BUILD)gcm_io_h.o			\
	$(BUILD)gefs_io.o          	\
	$(BUILD)gefs_io_h.o			\
	$(BUILD)obs_io.o    		\
	$(BUILD)obs_io_h.o			\
	$(BUILD)geo_reader.o      	\
    $(BUILD)basic_stats.o     	\
	$(BUILD)random.o            \
	$(BUILD)sampling.o			\
	$(BUILD)sampling_h.o		\
	$(BUILD)sorting.o			\
	$(BUILD)quantile_mapping.o  \
	$(BUILD)analogs.o       	\
	$(BUILD)analogs_h.o       	\
	$(BUILD)downscaling.o       \
	$(BUILD)downscaling_h.o     \
	$(BUILD)regression.o		\
	$(BUILD)regression_h.o		\
	$(BUILD)output.o    		\
	$(BUILD)output_h.o    		\
	$(BUILD)main.o


CONFIG_TEST_OBJS =    			\
    $(BUILD)configuration.o  	\
	$(BUILD)data_structures.o	\
	$(BUILD)model_constants.o	\
	$(BUILD)io_routines.o     	\
	$(BUILD)io_routines_h.o    	\
	$(BUILD)string.o     		\
	$(BUILD)time.o				\
	$(BUILD)test_config.o

CAL_TEST_OBJS =    			\
    $(BUILD)time.o 			\
	$(BUILD)test_calendar.o

SORT_TEST_OBJS = 			\
	$(BUILD)sorting.o		\
	$(BUILD)test_sort.o

QM_TEST_OBJS = 					\
	$(BUILD)data_structures.o	\
	$(BUILD)basic_stats.o		\
	$(BUILD)quantile_mapping.o	\
	$(BUILD)sorting.o			\
	$(BUILD)test_qm.o

REGRESSION_TEST_OBJS = 			\
	$(BUILD)regression.o		\
	$(BUILD)regression_h.o		\
	$(BUILD)test_regression.o

RANDOM_TEST_OBJS = 				\
	$(BUILD)test_random.o		\
	$(BUILD)random.o			\
	$(BUILD)model_constants.o

# These are the primary user facing rules
#
# make [all, test, install, clean, allclean]
#
# "make" builds model
model: gard

all: gard test

test: test_calendar test_sort test_qm test_config test_random  #test_regression

install: gard
	cp gard $(INSTALL_DIR)/

clean:
	$(RM) $(BUILD)* 2>/dev/null ||:

allclean:clean
	$(RM) gard test_calendar test_sort test_qm test_config test_regression 2>/dev/null ||:
	$(RM) -r *.dSYM 2>/dev/null ||:

# higher level rules
gard: $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

# various unit tests
# test_regression : $(REGRESSION_TEST_OBJS)
# 	$(FC) -o $@ $^ $(LDFLAGS)

test_calendar : $(CAL_TEST_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

test_qm : $(QM_TEST_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

test_sort : $(SORT_TEST_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

test_random : $(RANDOM_TEST_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

test_config : $(CONFIG_TEST_OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)



# specify specific builds and dependencies.

# unit tests

# The regression unit test is currently broken because it relied on a non-public API
# $(BUILD)test_regression.o: tests/test_regression.f90 $(BUILD)regression_h.o
# 	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)test_calendar.o: tests/test_calendar.f90 $(BUILD)time.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)test_sort.o: tests/test_sort.f90 $(BUILD)sorting.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)test_qm.o: tests/test_qm.f90 $(BUILD)data_structures.o $(BUILD)quantile_mapping.o $(BUILD)basic_stats.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)test_config.o: tests/test_config.f90 $(BUILD)data_structures.o $(BUILD)model_constants.o $(BUILD)configuration.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)test_random.o: tests/test_random.f90 $(BUILD)random.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@


# Date and time relevant objects
$(BUILD)time.o: timing/time.f90 $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)time_utils.o: timing/time_utils.f90 $(BUILD)time.o $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)time_io.o: timing/time_io.f90 $(BUILD)time.o $(BUILD)io_routines_h.o \
					$(BUILD)model_constants.o $(BUILD)string.o $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

# other utilities
$(BUILD)quantile_mapping.o: utilities/quantile_map.f90 $(BUILD)data_structures.o $(BUILD)sorting.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)sorting.o: utilities/sorting.f90 utilities/qsort_inline.inc
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)random.o: utilities/random.f90 $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)sampling_h.o: downscaling_stats/sampling_h.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)sampling.o: downscaling_stats/sampling.f90 $(BUILD)random.o $(BUILD)sampling_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)basic_stats.o: utilities/basic_stats.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)string.o: utilities/string.f90 $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

# constants and data structures
$(BUILD)model_constants.o: objects/model_constants.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)data_structures.o: objects/data_structures.f90 $(BUILD)time.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

# Input, output / configuration routines
$(BUILD)io_routines_h.o: io/io_routines_h.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)io_routines.o: io/io_routines.f90 $(BUILD)model_constants.o $(BUILD)io_routines_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)configuration.o: config/configuration.f90 $(BUILD)data_structures.o $(BUILD)model_constants.o \
					$(BUILD)string.o $(BUILD)io_routines_h.o $(BUILD)time.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)geo_reader.o: utilities/geo_reader.f90 $(BUILD)data_structures.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)init.o: main/init.f90 $(BUILD)data_structures.o $(BUILD)gcm_io_h.o $(BUILD)gefs_io_h.o \
					$(BUILD)obs_io_h.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)gefs_io_h.o: io/gefs_io_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)gefs_io.o: io/gefs_io.f90 $(BUILD)gefs_io_h.o $(BUILD)string.o $(BUILD)model_constants.o \
						$(BUILD)time_io.o $(BUILD)basic_stats.o $(BUILD)geo_reader.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)gcm_io_h.o: io/gcm_io_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)gcm_io.o: io/gcm_io.f90 $(BUILD)gcm_io_h.o $(BUILD)string.o $(BUILD)model_constants.o \
					$(BUILD)time_io.o $(BUILD)basic_stats.o $(BUILD)geo_reader.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)obs_io_h.o: io/obs_io_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)obs_io.o: io/obs_io.f90 $(BUILD)obs_io_h.o $(BUILD)string.o $(BUILD)model_constants.o \
					$(BUILD)time_io.o $(BUILD)basic_stats.o $(BUILD)geo_reader.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)output_h.o: io/output_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)output.o: io/output.f90 $(BUILD)output_h.o $(BUILD)io_routines_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)regression_h.o: downscaling_stats/regression_h.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

# Statical downscaling routines
$(BUILD)regression.o: downscaling_stats/regression.f90 $(BUILD)regression_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)atmosphere_io_h.o: io/atmosphere_io_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)atmosphere_io.o: io/atmosphere_io.f90 $(BUILD)atmosphere_io_h.o  \
					$(BUILD)gcm_io_h.o $(BUILD)gefs_io_h.o $(BUILD)model_constants.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)downscaling.o: downscaling_stats/downscaling.f90 $(BUILD)model_constants.o $(BUILD)downscaling_h.o \
					$(BUILD)quantile_mapping.o $(BUILD)regression_h.o $(BUILD)analogs_h.o $(BUILD)random.o $(BUILD)sampling_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)downscaling_h.o: downscaling_stats/downscaling_h.f90 $(BUILD)data_structures.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)analogs_h.o: downscaling_stats/analogs_h.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)analogs.o: downscaling_stats/analogs.f90 $(BUILD)string.o $(BUILD)analogs_h.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@

$(BUILD)main.o: main/main.f90 	$(BUILD)obs_io.o    								\
					$(BUILD)atmosphere_io_h.o $(BUILD)configuration.o				\
					$(BUILD)init.o $(BUILD)model_constants.o $(BUILD)io_routines_h.o\
					$(BUILD)downscaling_h.o $(BUILD)output_h.o $(BUILD)geo_reader.o
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@


$(BUILD)%.o: %.f90
	$(FC) $(FCFLAGS) -c $< $(MODOUTPUT) -o $@
