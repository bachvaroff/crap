# CC=/bin/time /bin/cc
CC=/bin/cc
CFLAGS=-D__BRANCH
CDEBUG=-D_DEBUG_ -D_DEBUG_CI_ -D_DEBUG_CIJ_
INCLUDE=-I.
LDFLAGS=
LDEBUG=
LIB=
LDPOSTFLAGS0=-lm
SOFLAGS=
SODEBUG=
COMPILE=$(CC) $(INCLUDE) $(CFLAGS) $(CDEBUG)
LINK=$(CC) $(LIB) $(LDFLAGS) $(LDEBUG)
LINKSO=$(CC) $(LIB) $(SOFLAGS) $(SODEBUG)

.c.o:
	$(COMPILE) -c $< -o $@

all: mult rot ttrig thyp tinner

cmatrix.o: cmatrix.c complex.h cmatrix.h
mult.o: mult.c complex.h cmatrix.h
rot.o: rot.c complex.h
ttrig.o: ttrig.c complex.h
thyp.o: thyp.c complex.h
tinner.o: tinner.c complex.h

mult: mult.o cmatrix.o
	$(LINK) mult.o cmatrix.o $(LDPOSTFLAGS0) -o $@

rot: rot.o
	$(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o
	$(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

thyp: thyp.o
	$(LINK) thyp.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o
	$(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o mult rot ttrig thyp tinner
