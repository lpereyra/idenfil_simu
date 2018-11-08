#CC
CC     := cbuild.csh

.PHONY : clean

MAKEFILE := Makefile

HEADERS := $(patsubst %.o,$.h,$(OBJS))

EXEC := pos

todo: $(EXEC)

pos: pos.c
	$(CC) $^ 

clean:
	rm -rf $(EXEC)
