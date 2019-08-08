### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE         #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in simulation Mpc/h
#EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
#EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DMCRITIC
#EXTRAS += -DWRITE_WEIGHT

#CC
CC     := $(OMPP) g++ $(DOMPP)
DC     := -DNTHREADS=32
#DC     += -DLOCK
CFLAGS := -Wall -O3 -fopenmp -g
#GSLL   := -lgsl -lgslcblas
VPP_INC := -I/home/lpereyra/voro++/include/voro++
VPP_LIB := -L/home/lpereyra/voro++/lib -lvoro++
#VPP_INC := -I/usr/local/include/voro++
#VPP_LIB := -lvoro++
LIBS   := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o leesnap.o grid.o voronoi.o kruskal.o

HEADERS := $(patsubst %.o,$.hh,$(OBJS))

EXEC := mst.x

todo: $(EXEC)

%.o: %.cc %.hh $(MAKEFILE)
	$(CC) $(EXTRAS) $(VPP_INC) $(CFLAGS) $(DC) -c $<

mst.x: mst.cc $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(VPP_LIB) $(CFLAGS) $(DC) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf mst.o
	rm -rf $(EXEC)
