### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
#EXTRAS += -DMPC                #Positions in Mpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DSTORE_IDS
EXTRAS += -DHUECOS
#EXTRAS += -DMCRITIC

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=16
#DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o grid.o variables.o calcula_random.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := cilindros_random.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

cilindros_random.x: cilindros_random.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf cilindros_random.o
	rm -rf $(EXEC)
