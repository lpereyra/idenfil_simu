### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DSTORE_IDS
EXTRAS += -DCALCULA_MEDIA
EXTRAS += -DHUECOS
EXTRAS += -DEXTEND
#EXTRAS += -DFIXED_SEPARATION
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DBIN_LOG
#EXTRAS += -DSAVEPART
#EXTRAS += -DMCRITIC
EXTRAS += -DREORDER
#EXTRAS += -DCOLUMN

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=16
CFLAGS := -Wall -O3 -fopenmp -march=native -ftree-vectorize -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o grid.o variables.o calcula.o list.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := cilindros.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

cilindros.x: cilindros.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf cilindros.o
	rm -rf $(EXEC)
