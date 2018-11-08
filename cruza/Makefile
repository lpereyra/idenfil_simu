### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
#EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DSTORE_IDS

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=56
CFLAGS := -Wall -O3 -fopenmp -march=native -g
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := cruza.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

cruza.x: cruza.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf cruza.o
	rm -rf $(EXEC)
