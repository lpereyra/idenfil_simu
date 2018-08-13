### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
#EXTRAS += -DEXTEND
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DBIN_LOG
#EXTRAS += -DMCRITIC
#EXTRAS += -DCALCULA_MEDIA
#EXTRAS += -DFIX_NSMOOTH
EXTRAS += -DITERA
#EXTRAS += -DNEW

#CC
CC      := $(OMPP) gcc $(DOMPP)
DC      := -DNTHREADS=4
CFLAGS  := -Wall -O3 -fopenmp -march=native -g
LDFLAGS := -lm -lgfortran 
FC      := gfortran
FFLAGS  := -cpp -O3 -fbounds-check 

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o libfitpack.o

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

EXEC := smooth.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

%.f.o: $(OBJS_FITPACK) 
	$(FC) $(FFLAGS) -c $<

smooth.x: smooth.c $(OBJS) $(OBJS_FITPACK)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJS_FITPACK)
	rm -rf $(OBJS)
	rm -rf smooth.o
	rm -rf $(EXEC)
