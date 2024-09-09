CC=cc
CFLAGS=-O -M3 -D_USE_ASM_ -D_USE_BRANCH_
CDEBUG=
INCLUDE=-I.
AS=cc
LD=cc
LDFLAGS=-M3
LDEBUG=
AR=ar
LIB=libcomplex.a
LIBH=cconfig.h complex.h cmatrix.h clu.h cfft.h fasm.h
LIBO=cmatrix.o clu.o cfft.o fasm.o
LDPOSTFLAGS0=$(LIB) -lm
SOFLAGS=
SODEBUG=
COMPILE=$(CFLAGS) $(CDEBUG) $(INCLUDE)
LINK=$(LDFLAGS) $(LDEBUG)
LINKSO=$(SOFLAGS) $(SODEBUG)

.c.o:
	$(CC) -c $(COMPILE) $< -o $@

.s.o:
	$(AS) -c $(COMPILE) $< -o $@

all: lib mult rot ttrig thyp tinner tlu inv tfft

cmatrix.o: cmatrix.c $(LIBH)
clu.o: clu.c $(LIBH)
cfft.o: cfft.c $(LIBH)
$(LIB): $(LIBO)
	$(AR) ruv $@ $?
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
	$(LD) $(LINK) mult.o $(LDPOSTFLAGS0) -o $@

inv: inv.o $(LIB)
	$(LD) $(LINK) inv.o $(LDPOSTFLAGS0) -o $@

rot: rot.o $(LIB)
	$(LD) $(LINK) rot.o $(LDPOSTFLAGS0) -o $@

ttrig: ttrig.o $(LIB)
	$(LD) $(LINK) ttrig.o $(LDPOSTFLAGS0) -o $@

thyp: thyp.o $(LIB)
	$(LD) $(LINK) thyp.o $(LDPOSTFLAGS0) -o $@

tinner: tinner.o $(LIB)
	$(LD) $(LINK) tinner.o $(LDPOSTFLAGS0) -o $@

tlu: tlu.o $(LIB)
	$(LD) $(LINK) tlu.o $(LDPOSTFLAGS0) -o $@

tfft: tfft.o $(LIB)
	$(LD) $(LINK) tfft.o $(LDPOSTFLAGS0) -o $@

clean:
	rm -f *.o $(LIB) mult rot ttrig thyp tinner tlu inv tfft

