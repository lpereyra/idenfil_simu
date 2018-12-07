### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE         #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in simulation Mpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
#EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DMCRITIC

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=30
#DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o leesnap.o chull.o

HEADERS := $(patsubst %.o,$.hh,$(OBJS))

EXEC := convex_hull.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

convex_hull.x: convex_hull.c $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(CFLAGS) $(DC) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf convex_hull.o
	rm -rf $(EXEC)
