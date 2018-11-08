### snapshot options #######
#EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
#EXTRAS += -DSTORE_IDS
#EXTRAS += -DSAVENEIGH
#EXTRAS += -DPOSFACTOR=1.
EXTRAS += -DSAVECENTROID
EXTRAS += -DWRITE_ASCII
EXTRAS += -DLLOYD
EXTRAS += -DNX=4
EXTRAS += -DNY=4
EXTRAS += -DNZ=2
EXTRAS += -DNTHREADS=32

#CC
CC     := g++
CFLAGS := -Wall -O3 -march=native -g
OMPP   := -fopenmp
VPP_DIR := /home/lpereyra/voro++
VPP_INC := -I$(VPP_DIR)/include/voro++
VPP_LIB := -L$(VPP_DIR)/lib -lvoro++
#VPP_DIR := /usr/local
#VPP_INC := -I$(VPP_DIR)/include/voro++
#VPP_LIB := -lvoro++
LIBS   := -lm 

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o leesnap.o voronoi.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := main.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(VPP_INC) $(CFLAGS) $(OMPP) -c $<

main.x: main.c $(OBJS)
	$(CC) $(EXTRAS) $^ $(LIBS) $(VPP_LIB) $(CFLAGS) $(OMPP) -o $@
	
clean:
	rm -rf $(OBJS)
	rm -rf main.o
	rm -rf $(EXEC)
