### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS           #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
#EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DBIN_LOG
#EXTRAS += -DCUT_ELONGACION
EXTRAS += -DCUT_IN_LEN
EXTRAS += -DTYPE_FLAG=1
#EXTRAS += -DTYPE_FLAG=2
EXTRAS += -DSIN_REPETICION
#EXTRAS += -DSERIAL
#EXTRAS += -DESFERA
#EXTRAS += -DMARIO
#EXTRAS += -DSINTETICOS_RECTOS

#CC
CC      := $(OMPP) mpicc $(DOMPP)
CFLAGS  := -Wall -O2 -fopenmp -march=native -ftree-vectorize -msse2 -ftree-vectorizer-verbose=5 -g
LDFLAGS := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o grid.o calcula.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := density_profile.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) -c $<

%.f.o: $(OBJS_FITPACK) 
	$(FC) $(FFLAGS) -c $<

density_profile.x: density_profile.c $(OBJS) $(OBJS_FITPACK)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJS)
	rm -rf density_profile.o
	rm -rf $(EXEC)
