# installation for amberlite:

include ../config.h

# trick to use $(AMBER_PREFIX)/bin
export AMBERHOME=$(AMBER_PREFIX)

install: $(BINDIR)/nab$(SFX) $(BINDIR)/ffgbsa$(SFX) $(BINDIR)/minab$(SFX) \
      $(BINDIR)/mdnab$(SFX)
	(if [ "$(SKIP_PYTHON)" = "no" ]; then \
		$(PYTHON) setup.py install -f $(PYTHON_INSTALL) ;\
	fi;\
	)

$(BINDIR)/nab$(SFX): ../nab/nab.c
	cd ../nab && $(MAKE) $@

$(BINDIR)/ffgbsa$(SFX): ffgbsa.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/ffgbsa$(SFX) ffgbsa.nab

$(BINDIR)/minab$(SFX): minab.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/minab$(SFX) minab.nab

$(BINDIR)/mdnab$(SFX): mdnab.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/mdnab$(SFX) mdnab.nab

clean:
	/bin/rm -f *.c

uninstall:
	(cd $(BINDIR) && /bin/rm -f ffgbsa minab mdnab pytleap pymdpbsa pdb4amber)
