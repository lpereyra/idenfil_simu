### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DMCRITIC
#EXTRAS += -DORIGINAL
#EXTRAS += -DPARTICLES

#CXX
CXX     := $(OMPP) g++ $(DOMPP)
DC     := -DNTHREADS=32
DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -march=native -g
#GSLL   := -lgsl -lgslcblas
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o grid.o list.o calcula.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := esferas_halos_filaments.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CXX) $(EXTRAS) $(CFLAGS) $(DC) -c $<

esferas_halos_filaments.x: esferas_halos_filaments.c $(OBJS)
	$(CXX) $(EXTRAS) $^ $(CFLAGS) $(DC) $(LIBS) -o $@ 

clean:
	rm -rf $(OBJS)
	rm -rf esferas_halos_filaments.o
	rm -rf $(EXEC)
