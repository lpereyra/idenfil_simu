### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DCUT_IN_LEN
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DMCRITIC
#EXTRAS += -DSPH 
#EXTRAS += -DRAUX=1000. 
#EXTRAS += -DCENTERS_RANDOM

#CXX
CXX     := $(OMPP) g++ $(DOMPP)
DC     := -DNTHREADS=32
DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -march=native -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o calcula_pares_random.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := new_pares_random.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CXX) $(EXTRAS) $(CFLAGS) $(DC) -c $<

new_pares_random.x: new_pares_random.c $(OBJS)
	$(CXX) $(EXTRAS) $^ $(CFLAGS) $(DC) $(LIBS) -o $@ 

clean:
	rm -rf $(OBJS)
	rm -rf new_pares_random.o
	rm -rf $(EXEC)
