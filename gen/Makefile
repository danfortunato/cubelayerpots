MWRAP     := mwrap
MEX       := mex
FC        := gfortran
FFLAGS    := -O3 -march=native -fopenmp -std=legacy -w -fPIC
MEXFLAGS  := -DMWF77_UNDERSCORE1 LDFLAGS='$$LDFLAGS -mmacosx-version-min=15.0'
MEXLIBS   := -L/opt/homebrew/lib/gcc/current -lgfortran -lgomp

MEX_NAME  := gentab_mex
MW_FILE   := gentab.mw
M_FILE    := gentab.m
C_GATEWAY := $(MEX_NAME).c
F_SRCS    := gentab_omp.f90 \
           	 genslpsame_omp.f90 \
           	 genslpside_omp.f90 \
           	 gendlpside_omp.f90 \
           	 dcuhre_nocommon.f \
           	 chebexps.f \
           	 prini.f
F_OBJS    := $(patsubst %.f90,%.o,$(patsubst %.f,%.o,$(F_SRCS)))

.PHONY: all mwrap mex test clean

all: mex

mwrap: $(C_GATEWAY) $(M_FILE)

$(C_GATEWAY) $(M_FILE): $(MW_FILE)
	$(MWRAP) -mex $(MEX_NAME) -m $(M_FILE) -c $(C_GATEWAY) $(MW_FILE)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

mex: $(MEX_NAME).mex*

$(MEX_NAME).mex*: $(C_GATEWAY) $(F_OBJS)
	$(MEX) $(MEXFLAGS) $(C_GATEWAY) $(F_OBJS) $(MEXLIBS)

clean:
	rm -f $(F_OBJS)
	rm -f $(MEX_NAME).mex*

mwrapclean: clean
	rm -f $(C_GATEWAY)
	rm -f $(M_FILE)
