CC=cc
CFLAGS=-O -M3 -D_USE_ASM_ -D_USE_BRANCH_
CDEBUG=
INCLUDE=-I.
LDFLAGS=-M3
LDEBUG=
LIB=libcomplex.a
LIBH=complex.h cmatrix.h clu.h cfft.h
LDPOSTFLAGS0=$(LIB) -lm
SOFLAGS=
SODEBUG=
COMPILE=$(CFLAGS) $(CDEBUG) $(INCLUDE)
LINK=$(LDFLAGS) $(LDEBUG)
LINKSO=$(SOFLAGS) $(SODEBUG)

.c.o:
	$(CC) -c $(COMPILE) $< -o $@

.s.o:
	$(CC) -c $(COMPILE) $< -o $@

all: lib mult rot ttrig thyp tinner tlu inv tfft

cmatrix.o: cmatrix.c cconfig.h complex.h fasm.h cmatrix.h
clu.o: clu.c cconfig.h complex.h fasm.h clu.h
cfft.o: cfft.c cconfig.h complex.h fasm.h cfft.h
$(LIB): cmatrix.o clu.o cfft.o fasm.o
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
	$(CC) $(LINK) mult.o $(LDPOSTFLAGS0) -o $@

inv: inv.o $(LIB)
	$(CC) $(LINK) inv.o $(LDPOSTFLAGS0) -o $@

rot: rot.o $(LIB)
	$(CC) $(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o $(LIB)
	$(CC) $(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

thyp: thyp.o $(LIB)
	$(CC) $(LINK) thyp.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o $(LIB)
	$(CC) $(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

tlu: tlu.o $(LIB)
	$(CC) $(LINK) tlu.o $(LDPOSTFLAGS0) -o $@

tfft: tfft.o $(LIB)
	$(CC) $(LINK) tfft.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o $(LIB) mult rot ttrig thyp tinner tlu inv tfft
