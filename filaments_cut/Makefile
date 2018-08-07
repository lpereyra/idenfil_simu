### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE         #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
#EXTRAS += -DPOSFACTOR=1000.0   #Positions in simulation Mpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DPRUNED
EXTRAS += -DBRANCH_SURVIVE
EXTRAS += -DSORT_DERECHA
#EXTRAS += -DLEVEL_PRUNED=10
#EXTRAS += -DMCRITIC
EXTRAS += -DNEW

#CC
CC     := $(OMPP) g++ $(DOMPP)
DC     := -DNTHREADS=16
CFLAGS := -Wall -O3 -fopenmp -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o leesnap.o pruned.o

HEADERS := $(patsubst %.o,$.hh,$(OBJS))

EXEC := fil.x

todo: $(EXEC)

%.o: %.cc %.hh $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

fil.x: fil.cc $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(CFLAGS) $(DC) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf fil.o
	rm -rf $(EXEC)
