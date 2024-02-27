CC=cc
CFLAGS=
INCLUDE=-I.
LDFLAGS=
LDPOSTFLAGS0=-lm
SOFLAGS=
DEBUG=-D_DEBUG_ -D_DEBUG_CI_ -D_DEBUG_CIJ_
COMPILE=$(CC) $(DEBUG) $(INCLUDE) $(CFLAGS)
LINK=$(CC) $(DEBUG) $(LDFLAGS)
LINKSO=$(CC) $(DEBUG) $(SOFLAGS)
LD_LIBRARY_PATH=

.c.o:
	$(COMPILE) -c $< -o $@

all: mult rot

mult: mult.o cmatrix.o
	$(LINK) mult.o cmatrix.o $(LDPOSTFLAGS0) -o $@

rot: rot.o
	$(LINK) rot.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o mult rot
