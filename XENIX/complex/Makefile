CC=cc
CFLAGS=-D__BRANCH
CDEBUG=-g -D_DEBUG_ -D_DEBUG_CI_
INCLUDE=-I.
LDFLAGS=
LDEBUG=-g
LIB=libcomplex.a
LIBH=complex.h cmatrix.h clu.h cfft.h
LDPOSTFLAGS0=$(LIB) -lm
SOFLAGS=
SODEBUG=-g
COMPILE=$(CC) $(CFLAGS) $(CDEBUG) $(INCLUDE)
LINK=$(CC) $(LDFLAGS) $(LDEBUG)
LINKSO=$(CC) $(SOFLAGS) $(SODEBUG)

.c.o:
	$(COMPILE) -c $< -o $@

all: lib mult rot ttrig thyp tinner tlu inv tfft

cmatrix.o: cmatrix.c complex.h cmatrix.h
clu.o: clu.c complex.h clu.h
cfft.o: cfft.c complex.h cfft.h
$(LIB): cmatrix.o clu.o cfft.o
	ar ruv $@ $?
	ranlib $@
lib: $(LIB)

mult.o: mult.c $(LIBH)
rot.o: rot.c $(LIBH)
ttrig.o: ttrig.c $(LIBH)
thyp.o: thyp.c $(LIBH)
tinner.o: tinner.c $(LIBH)
tlu.o: tlu.c $(LIBH)
inv.o: inv.c $(LIBH)
tfft.o: tfft.c $(LIBH)

mult: mult.o $(LIB)
	$(LINK) mult.o $(LDPOSTFLAGS0) -o $@

inv: inv.o $(LIB)
	$(LINK) inv.o $(LDPOSTFLAGS0) -o $@

rot: rot.o $(LIB)
	$(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o $(LIB)
	$(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

thyp: thyp.o $(LIB)
	$(LINK) thyp.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o $(LIB)
	$(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

tlu: tlu.o $(LIB)
	$(LINK) tlu.o $(LDPOSTFLAGS0) -o $@

tfft: tfft.o $(LIB)
	$(LINK) tfft.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o $(LIB) mult rot ttrig thyp tinner tlu inv tfft
