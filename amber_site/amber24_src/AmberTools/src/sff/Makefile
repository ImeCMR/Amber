include ../config.h

.c.o:
	@echo "[SFF]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERCFLAGS) $(RISMSFF) $(NETCDFINC) -o $@ $<

OBJS = conjgrad.o lmodC.o memutil.o nblist.o newton.o nmode.o \
	   prm.o rand2.o sasad.o sff.o xminC.o AmberNetcdf.o ptrajmask.o \
       $(SFF_RISM_INTERFACE)

install: libsff hcp_getpdb$(SFX)

libsff: $(OBJS)
	@echo "[SFF]  AR $@"
	$(VB)/bin/rm -f $(LIBDIR)/libsff.a
	$(VB)$(AR) $(LIBDIR)/libsff.a $(OBJS)	> /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libsff.a
	$(VB)cp sff.h AmberNetcdf.h $(INCDIR)

libsff_mpi: $(OBJS)
	@echo "[SFF]  AR $@"
	$(VB)/bin/rm -f $(LIBDIR)/libsff_mpi.a
	$(VB)$(AR) $(LIBDIR)/libsff_mpi.a $(OBJS)	> /dev/null
	$(VB)$(RANLIB) $(LIBDIR)/libsff_mpi.a
	$(VB)cp sff.h AmberNetcdf.h $(INCDIR)

uninstall:
	/bin/rm -f $(BINDIR)/hcp_getpdb$(SFX)
	/bin/rm -f $(LIBDIR)/libsff.a $(LIBDIR)/libsff_mpi.a
	/bin/rm -f $(INCDIR)/sff.h $(INCDIR)/AmberNetcdf.h

clean:
	/bin/rm -f $(OBJS) lex.mmo.c

lex.mmo.c: mm_options.l
	@echo "[SFF]  FLEX $<"
	$(VB)$(LEX) mm_options.l

hcp_getpdb$(SFX): hcp_getpdb.c
	@echo "[SFF]  CC $<"
	$(VB)$(CC) -o $(BINDIR)/hcp_getpdb$(SFX) hcp_getpdb.c

.PHONY: ../rism/amber_rism_interface.NAB.o
../rism/amber_rism_interface.NAB.o:
	cd ../rism && \
	$(MAKE) amber_rism_interface.NAB.o

sff.c: lex.mmo.c

depend: conjgrad.c lmodC.c memutil.c nblist.c newton.c nmode.c \
       prm.c rand2.c sasad.c sff.c xminC.c AmberNetcdf.c ptrajmask.c
	$(CC) -MM $^ > $@

include depend

