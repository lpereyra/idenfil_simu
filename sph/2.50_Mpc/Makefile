### snapshot options #######
#EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
#EXTRAS += -DSTORE_IDS
#EXTRAS += -DSAVENEIGH
#EXTRAS += -DPOSFACTOR=1.
EXTRAS += -DCUT_IN_LEN
EXTRAS += -DLEN_MIN=4000.
EXTRAS += -DLEN_MAX=101000.
#EXTRAS += -DWRITE_ASCII

#CC
CC      := $(OMPP) g++ $(DOMPP)
DC      := -DNTHREADS=4
CFLAGS  := -Wall -O3 -fopenmp -march=native -g
LDFLAGS := -lm -lgfortran 
FF      := gfortran
FFLAGS  := -cpp -O3 -fbounds-check -g -std=legacy

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := variables.o grid.o leesnap.o sph.o main.o 
OBJS_ARVO := arvo.o 
						  
HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := main.x

todo: $(EXEC)

%.o: %.f
	$(FF) $(FFLAGS) -c $<

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

main.x: $(OBJS_ARVO) $(OBJS)
	$(CC) $(EXTRAS) $(CFLAGS) $(OMPP) $(DC) $^ -o $(EXEC) $(LDFLAGS) 

clean:
	rm -rf $(OBJS_ARVO) $(OBJS)
	rm -rf main.o
	rm -rf $(EXEC)
