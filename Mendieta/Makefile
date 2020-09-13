### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0    #Positions in Kpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
#EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DSTORE_IDS
#EXTRAS += -DCOLUMN
EXTRAS += -DFILE_ASCII
#EXTRAS += -DCHANGE_POSITION

#CC
CC     := $(OMPP) gcc $(DOMPP)
DC     := -DNTHREADS=4
DC     += -DLOCK
GSLL   := -lgsl -lgslcblas
CFLAGS := -Wall -O3 -march=native -ftree-vectorize -fopenmp -g
LIBS   := -lm $(GSLL) 

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS :=  variables.o io.o allocate.o leesnap.o grid.o iden.o propiedades.o deltas.o octree.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mendieta.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

mendieta.x: mendieta.c $(OBJS)
	$(CC) $(CFLAGS) $(EXTRAS) $^  -o $@ $(LIBS)

clean:
	rm -rf $(OBJS)
	rm -rf mendieta.o
	rm -rf $(EXEC)
