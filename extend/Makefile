### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DORIGINAL
#EXTRAS += -DNEW
EXTRAS += -DLAST_VECTOR

#CC
CC      := $(OMPP) gcc $(DOMPP)
DC      := -DNTHREADS=8
CFLAGS  := -Wall -O3 -fopenmp -march=native -g
LDFLAGS := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o 

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := extend.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

extend.x: extend.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJS)
	rm -rf extend.o
	rm -rf $(EXEC)
