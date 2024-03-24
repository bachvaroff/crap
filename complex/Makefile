CC=/bin/cc
CFLAGS=-D__BRANCH
CDEBUG=-g -D_DEBUG_ -D_DEBUG_CI_
# CDEBUG=-D_DEBUG_ -D_DEBUG_CI_ -D_DEBUG_CIJ_
INCLUDE=-I.
LDFLAGS=
LDEBUG=-g
LIB=
LDPOSTFLAGS0=-lm
SOFLAGS=
SODEBUG=-g
COMPILE=$(CC) $(CFLAGS) $(CDEBUG) $(INCLUDE)
LINK=$(CC) $(LDFLAGS) $(LDEBUG) $(LIB)
LINKSO=$(CC) $(SOFLAGS) $(SODEBUG) $(LIB)

CMOBJ=cmatrix.o clu.o

.c.o:
	$(COMPILE) -c $< -o $@

all: mult rot ttrig thyp tinner tlu inv

cmatrix.o: cmatrix.c complex.h cmatrix.h
clu.o: clu.c complex.h cmatrix.h clu.h
mult.o: mult.c complex.h cmatrix.h
rot.o: rot.c complex.h
ttrig.o: ttrig.c complex.h
thyp.o: thyp.c complex.h
tinner.o: tinner.c complex.h
tlu.o: tlu.c complex.h cmatrix.h clu.h
inv.o: inv.c complex.h cmatrix.h clu.h

mult: mult.o $(CMOBJ)
	$(LINK) mult.o $(CMOBJ) $(LDPOSTFLAGS0) -o $@

inv: inv.o $(CMOBJ)
	$(LINK) inv.o $(CMOBJ) $(LDPOSTFLAGS0) -o $@

rot: rot.o
	$(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o
	$(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

thyp: thyp.o
	$(LINK) thyp.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o
	$(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

tlu: tlu.o $(CMOBJ)
	$(LINK) tlu.o $(CMOBJ) $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o mult rot ttrig thyp tinner tlu inv
