### snapshot options #######
EXTRAS += -DPERIODIC    #periodic boundary condition
#EXTRAS += -DPRECDOUBLE   #Pos and vel in double precision
#EXTRAS += -DLONGIDS            #IDs are long integer
EXTRAS += -DPOSFACTOR=1000.0   #Positions in Kpc
EXTRAS += -DVELFACTOR=1.0      #Velocities in km/s
EXTRAS += -DSTORE_VELOCITIES
#EXTRAS += -DCUT_IN_LEN
EXTRAS += -DHUECOS
EXTRAS += -DEXTEND
#EXTRAS += -DVEL_RELATIVA
#EXTRAS += -DBIN_LOG
#EXTRAS += -DMCRITIC
EXTRAS += -DCALCULA_VCM
#EXTRAS += -DCALCULA_MEDIA
#EXTRAS += -DNEW
#EXTRAS += -DORIGINAL
#EXTRAS += -DCUT_ELONGACION
EXTRAS += -DTYPE_FLAG=1

#CC
CC      := $(OMPP) gcc $(DOMPP)
DC      := -DNTHREADS=56
DC      += -DLOCK
CFLAGS  := -Wall -O3 -fopenmp -march=native -g
LDFLAGS := -lm

.PHONY : cleanall clean todo 

MAKEFILE := Makefile

OBJS := leesnap.o variables.o grid.o list.o calcula.o

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := new_cilindros_extend.x

todo: $(EXEC)

%.o: %.c %.h $(MAKEFILE)
	$(CC) $(EXTRAS) $(CFLAGS) $(DC) -c $<

%.f.o: $(OBJS_FITPACK) 
	$(FC) $(FFLAGS) -c $<

new_cilindros_extend.x: new_cilindros_extend.c $(OBJS) $(OBJS_FITPACK)
	$(CC) $(CFLAGS) $(EXTRAS) $(DC) $^  -o $@ $(LDFLAGS)

clean:
	rm -rf $(OBJS)
	rm -rf new_cilindros_extend.o
	rm -rf $(EXEC)
