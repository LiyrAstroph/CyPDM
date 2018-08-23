SHELL=/bin/bash
CC       = gcc 
OPTIMIZE = -O2 -Wall -finline-functions

EXEC     = pdm
SRC      = ./

OBJS = $(SRC)/pdm.o

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) 
LIBS     = -lm

INCL     = Makefile $(SRC)/pdm.h

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS) $(LIBS) -o $@

$(OBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC)