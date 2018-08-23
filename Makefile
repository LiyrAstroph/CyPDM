SHELL=/bin/bash
CC       = gcc 
OPTIMIZE = -O2 -Wall -finline-functions

LIBPDM  = libpdm.so
EXEC     = pdm
SRC      = ./

all: $(LIBPDM) $(EXEC)

LIBOBJS = $(SRC)/pdm.o 
OBJS = $(SRC)/main.o

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) 
LIBS     = -lm

INCL     = Makefile $(SRC)/pdm.h

$(LIBPDM): $(LIBOBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(LIBS) -fPIC -shared -o libpdm.so $(SRC)/pdm.c

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS) $(LIBS)  -L./ -lpdm -o $@

$(OBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC)
