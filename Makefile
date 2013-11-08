CC = gcc
AR = ar
OBJS = hull.o ch.o io.o rand.o pointops.o fg.o HDRS = hull.h points.h pointsites.h stormacs.h SRC = hull.c ch.c io.c rand.c pointops.c fg.c PROG = hull BINDIR = ./bin LIBDIR = ./lib # or just /usr/lib on a 32 bit machine?
LIB = $(LIBDIR)/lib$(PROG).a
 
all : $(PROG) rsites
	cp $(PROG) $(BINDIR)/.
	cp rsites $(BINDIR)/.
 
$(OBJS) : $(HDRS)
 
hullmain.o : $(HDRS)
 
$(PROG) : $(OBJS) hullmain.o
	$(CC) $(OBJS) hullmain.o -o $(PROG) -lm
	$(AR) rcv $(LIB) $(OBJS)
 
rsites : rsites.c
	$(CC) -o rsites rsites.c -lm
 
clean :
	-rm -f $(OBJS) hullmain.o core a.out $(PROG)
