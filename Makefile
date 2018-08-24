SHELL=/bin/bash
CC       = gcc 
OPTIMIZE = -O2 -Wall -finline-functions

LIBPDM  = libcpdm.so
EXEC     = cpdm
SRC      = ./

all: $(LIBPDM) $(EXEC)

LIBOBJS = $(SRC)/cpdm.o 
OBJS = $(SRC)/demo.o

OPTIONS  = $(OPTIMIZE)
CFLAGS   = $(OPTIONS) 
LIBS     = -lm

INCL     = Makefile $(SRC)/cpdm.h

$(LIBPDM): $(LIBOBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(LIBS) -fPIC -shared -o libcpdm.so $(SRC)/cpdm.c

$(EXEC): $(OBJS)
	cd $(SRC)
	$(CC) $(OPTIMIZE) $(CFLAGS) $(OBJS) $(LIBS)  -L./ -lcpdm -o $@

$(OBJS): $(INCL)

$(LIBOBJS): $(INCL)

clean:
	rm $(SRC)/*.o $(EXEC) $(LIBPDM)
