# Makefile for Python reimplementation of MM/PB(GB)SA written by Jason Swails
# (adapted a bit for installing to $AMBER_PREFIX by Hai)
include ../config.h
AMBERHOME=$(AMBER_PREFIX)

SHELL = /bin/sh

install: $(BINDIR)/mmpbsa_py_nabnmode$(SFX) $(BINDIR)/mmpbsa_py_energy$(SFX)
	$(PYTHON) setup.py install -f $(PYTHON_INSTALL) --install-scripts=$(BINDIR)

$(BINDIR)/mmpbsa_py_nabnmode$(SFX): mmpbsa_entropy.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/mmpbsa_py_nabnmode$(SFX) mmpbsa_entropy.nab

$(BINDIR)/mmpbsa_py_energy$(SFX): mmpbsa_energy.nab
	$(BINDIR)/nab$(SFX) -o $(BINDIR)/mmpbsa_py_energy$(SFX) mmpbsa_energy.nab

serial: install

parallel:
	-$(PYTHON) setup.py install -f $(PYTHON_INSTALL) --install-scripts=$(BINDIR)

clean:
	-(/bin/rm -fr build/ mmpbsa_entropy.c mmpbsa_energy.c)

uninstall:
	/bin/rm -fr $(BINDIR)/mmpbsa_py_nabnmode$(SFX) $(BINDIR)/mmpbsa_py_energy$(SFX)
	/bin/rm -fr $(BINDIR)/MMPBSA.py $(BINDIR)/MMPBSA.py.MPI $(BINDIR)/ante-MMPBSA.py
