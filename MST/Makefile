### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS           #IDs are long integer
EXTRAS += -DMPC                #Positions in simulation Mpc/h
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
EXTRAS += -DSTORE_IDS
EXTRAS += -DBRANCH_SURVIVE
EXTRAS += -DSORT_DERECHA
#EXTRAS += -DMCRITIC

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

OBJS := variables.o leesnap.o grid.o voronoi.o mst_kruskal.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := mst_new.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(VPP_INC) $(CFLAGS) $(DC) -c $<

mst_new.x: mst_new.c $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(VPP_LIB) $(CFLAGS) $(DC) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf mst_new.o
	rm -rf $(EXEC)
