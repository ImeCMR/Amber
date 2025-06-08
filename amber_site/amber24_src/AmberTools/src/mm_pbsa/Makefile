#************************************************************************
# Makefile for AmberTools24 MM_PBSA
#
SHELL=/bin/sh
include ../config.h

# trick to use $(AMBER_PREFIX)/bin/nab
export AMBERHOME=$(AMBER_PREFIX)

LIBSRC= ../lib/amopen.F ../lib/mexit.F

LIBOBJ= ../lib/amopen.o ../lib/mexit.o

PROGS = $(BINDIR)/make_crd_hg$(SFX) $(BINDIR)/mm_pbsa_nabnmode$(SFX)

all: $(PROGS)

$(BINDIR)/make_crd_hg$(SFX): make_crd_hg.o $(LIBOBJ)
	$(FC) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $(BINDIR)/make_crd_hg$(SFX) make_crd_hg.o \
		$(LIBOBJ) $(LDFLAGS) $(AMBERLDFLAGS)

#New nab-less mm_pbsa_nabnmode a work in progress
$(BINDIR)/mm_pbsa_nabnmode$(SFX): mm_pbsa_nabnmode.o

#Note dependency on nab for original mm_pbsa_nabnmode
$(BINDIR)/mm_pbsa_nabnmode.original$(SFX): mm_pbsa_nabnmode.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/mm_pbsa_nabnmode.original$(SFX) \
	$(NABFLAGS) mm_pbsa_nabnmode.nab

../lib/amopen.o: ../lib/amopen.F
	cd ../lib; $(MAKE) amopen.o

../lib/mexit.o: ../lib/mexit.F
	cd ../lib; $(MAKE) mexit.o

install: all
	cp mm_pbsa.pl mm_pbsa_statistics.pl $(BINDIR)
	cp *.pm $(LIBDIR)/perl/mm_pbsa/

clean:
	/bin/rm -f *.o

uninstall: clean
	/bin/rm -f $(BINDIR)/mm_pbsa.pl
	/bin/rm -f $(BINDIR)/mm_pbsa_statistics.pl
	-for file in $(PROGS) ; do \
		/bin/rm -f $${file} ; \
	done
