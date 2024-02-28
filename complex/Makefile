CC=/bin/time cc
CFLAGS=
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

cmatrix.o: complex.h cmatrix.h
mult.o: complex.h cmatrix.h
rot.o: complex.h
ttrig.o: complex.h
tinner.o: complex.h

all: mult rot ttrig tinner

mult: mult.o cmatrix.o
	$(LINK) mult.o cmatrix.o $(LDPOSTFLAGS0) -o $@

rot: rot.o
	$(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o
	$(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o
	$(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o mult rot ttrig tinner
