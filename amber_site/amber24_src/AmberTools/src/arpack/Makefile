SHELL=/bin/sh
SRCDIR=..
VPATH=
include ../config.h

.f.o:
	@echo "[ARPACK]  FC $<"
	$(VB)$(FC) -c $(FOPTFLAGS) $(FFLAGS) -o $@ $<

.c.o:
	@echo "[ARPACK]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) -o $@ $<

install: $(LIBDIR)/libarpack.a

ARPACKLIB= ivout.o dvout.o dgetv0.o dlaqrb.o dnaitr.o dnapps.o \
dnaup2.o dnaupd.o dnconv.o dneigh.o dmout.o \
dneupd.o dngets.o dsaitr.o dsapps.o dsaup2.o \
dsaupd.o dsconv.o dseigt.o dsesrt.o dseupd.o \
dsgets.o dsortc.o dsortr.o dstatn.o dstats.o dstqrb.o dsarpack.o \
wallclock.o

SARPACKLIB= svout.o sgetv0.o slaqrb.o snaitr.o snapps.o \
snaup2.o snaupd.o snconv.o sneigh.o smout.o \
sneupd.o sngets.o ssaitr.o ssapps.o ssaup2.o \
ssaupd.o ssconv.o sseigt.o ssesrt.o sseupd.o \
ssgets.o ssortc.o ssortr.o sstatn.o sstats.o sstqrb.o

$(LIBDIR)/libarpack.a: $(ARPACKLIB)
	@echo "[ARPACK]  CC $@"
	$(VB)$(AR) libarpack.a $(ARPACKLIB) > /dev/null
	$(VB)$(RANLIB) libarpack.a
	$(VB)mv libarpack.a $(LIBDIR)

clean:
	-/bin/rm -f *.o _*.f

uninstall:
	-rm -f $(LIBDIR)/libarpack.a

