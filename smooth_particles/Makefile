### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
#EXTRAS += -DEXTEND
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DBIN_LOG
#EXTRAS += -DMCRITIC
#EXTRAS += -DFIX_NSMOOTH
EXTRAS += -DITERA
EXTRAS += -DPARTICLES
#EXTRAS += -DCALCULA_VCM

#CC
CC      := $(OMPP) gcc $(DOMPP)
DC      := -DNTHREADS=56
CFLAGS  := -Wall -O3 -fopenmp -march=native -g
LDFLAGS := -lm -lgfortran 
FC      := gfortran
FFLAGS  := -cpp -O3 -fbounds-check 

.PHONY : clean todo 
 
MAKEFILE := Makefile

OBJS := leesnap.o variables.o libfitpack.o sph.o grid.o

OBJS_ARVO:= arvo.o

OBJS_FITPACK := ./fitpack/parcur.o \
							  ./fitpack/fpback.o \
						    ./fitpack/fpbspl.o \
						    ./fitpack/fpchec.o \
						    ./fitpack/fppara.o \
						    ./fitpack/fpdisc.o \
						    ./fitpack/fpgivs.o \
						    ./fitpack/fpknot.o \
						    ./fitpack/fprati.o \
						    ./fitpack/fprota.o \
							  ./fitpack/splev.o  

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := smooth_particles.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

%.f.o: $(OBJS_FITPACK) 
	$(FC) $(FFLAGS) -c $<

%.f.o: $(OBJS_ARVO) 
	$(FC) $(FFLAGS) -c $<


smooth_particles.x: smooth_particles.c $(OBJS) $(OBJS_FITPACK) $(OBJS_ARVO)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJS_FITPACK)
	rm -rf $(OBJS_ARVO)
	rm -rf $(OBJS)
	rm -rf smooth_particles.o
	rm -rf $(EXEC)
