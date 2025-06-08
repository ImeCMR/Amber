SHELL=/bin/sh
include ../config.h

install:
	$(CC) miscellaneous/draw_membrane2.c -lm -o miscellaneous/draw_membrane2
	AMBER_PREFIX=$(AMBER_PREFIX) AMBER_SOURCE=$(AMBER_SOURCE) /bin/sh setup.sh
	cp -f $(AMBER_SOURCE)/AmberTools/src/FEW/FEW.pl $(BINDIR)

clean:

uninstall:
	(cd $(BINDIR) && /bin/rm -f FEW.pl)
	/bin/rm -rf $(AMBERHOME)/AmberTools/src/FEW/additional_libs/PerlMol
