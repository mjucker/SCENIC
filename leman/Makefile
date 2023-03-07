FFTW=/usr/local/fftw
#FFTW=/usr/local/fftw-2.1.5-32

FFTWLIB=-L$(FFTW)/lib -lfftw
FFTWINC=-I$(FFTW)/include

PROG =	leman

SRCS =	erf.f fields.f inout.f math.f transition.f main.f matrix.f unicity.f walltime.f90 utils.c

OBJS =	erf.o fields.o inout.o math.o transition.o main.o matrix.o unicity.o walltime.o utils.o

LIBS =  -L/usr/local/intel/mkl/lib/em64t -lmkl_lapack \
	     -lmkl_em64t -lguide -lpthread

#LIBS =  -L/usr/local/intel/mkl72/lib/32 -lmkl_lapack \
#	     -lmkl_p4 -lguide -lpthread

CC = cc
CFLAGS = -O
FC = mpif90
FFLAGS = -w -O3 -xT -FR -vec_report1 $(FFTWINC) -openmp
#FFLAGS = -w -O3 -xW -FR -vec_report1 $(FFTWINC) 
#FFLAGS = -CB -CA -w -g -FR $(FFTWINC) -openmp
#FFLAGS = -g -traceback -CA -CB -O3 -xT -FR $(FFTWINC) -openmp
F90 = $(FC)
F90FLAGS = $(FFLAGS)
LDFLAGS =

all: $(PROG)

$(PROG): $(OBJS)
	$(FC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS) $(FFTWLIB)

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

math.o: cputime.inc
matrix.o: math.o
main.o: math.o matrix.o unicity.o utils.o
inout.o: math.o matrix.o
unicity.o: math.o matrix.o
fields.o: math.o matrix.o
