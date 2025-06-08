include ../config.h

ndiff$(SFX):	ndiff.o store.o awklib.o regexp.o regerror.o
	$(CC) $(CFLAGS) -o $@ ndiff.o store.o awklib.o regexp.o regerror.o

install:	ndiff$(SFX)
	/bin/mv ndiff$(SFX) $(BINDIR)

clean::
	/bin/rm -f *.o ndiff$(SFX)

uninstall:
	-rm -f $(BINDIR)/ndiff$(SFX)

ndiff.o:	ndiff.c awklib.h ndiff.h
awklib.o:	awklib.c awklib.h
store.o:	store.c ndiff.h
