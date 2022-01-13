goal:  EVEREST_sfit.exe

tarball:
	tar cf duo.tar makefile *.f90
        
checkin:
	ci -l makefile *.f90


FPATH = 

PLAT = 0101_ch

FOR  = ifort
FFLAGS =  -O3 -qopenmp -traceback  -ipo -mkl=parallel  -cpp 
#FFLAGS =  -O0 -C -check bounds -g  -gen-interfaces -warn interfaces  -check arg_temp_created -prof-value-profiling=all -warn all  -traceback  -mkl  -cpp -C 

#FOR = gfortran
#FFLAGS = -Wall -march=native -g -O3 -fopenmp -std=f2008 -ffree-line-length-512 -fall-intrinsics

###############################################################################

OBJ =  accuracy.o  timer.o input.o  xyz_pes_koput_cm.o

EVEREST_sfit.exe:	$(OBJ) EVEREST_sfit.o 
	$(FOR) -o j-ev-fit_$(PLAT).x $(OBJ) EVEREST_sfit.o  $(FFLAGS) $(LIB) -static

EVEREST_sfit.o:	EVEREST_sfit.f90 $(OBJ) 
	$(FOR) -c EVEREST_sfit.f90 $(FFLAGS)

xyz_pes_koput_cm.o:	xyz_pes_koput_cm.f90
	$(FOR) -c xyz_pes_koput_cm.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

clean:
	rm $(OBJ) *.mod






