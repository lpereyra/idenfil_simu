### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DMPC                #Positions in Mpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DSTORE_IDS
#EXTRAS += -DCALCULA_MEDIA
EXTRAS += -DHUECOS
#EXTRAS += -DFIXED_SEPARATION
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DBIN_LOG
#EXTRAS += -DSAVEPART
#EXTRAS += -DMCRITIC

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=16
#DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -march=native -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o calcula.o grid.o list.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := halos_filament.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

halos_filament.x: halos_filament.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf halos_filament.o
	rm -rf $(EXEC)
