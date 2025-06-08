include ../config.h

install:
	$(PYTHON) setup.py install --prefix=$(AMBERHOME) --install-scripts=$(BINDIR)

clean:
	/bin/rm -rf build/

uninstall:
	/bin/rm -f $(BINDIR)/MCPB.py $(BINDIR)/OptC4.py $(BINDIR)/PdbSearcher.py \
           $(BINDIR)/espgen.py $(BINDIR)/CartHess2FC.py $(BINDIR)/IPMach.py \
           $(BINDIR)/car_to_files.py $(BINDIR)/ProScrs.py $(BINDIR)/mol2rtf.py \
           $(BINDIR)/amb2chm_psf_crd.py $(BINDIR)/amb2chm_par.py $(BINDIR)/amb2gro_top_gro.py \
           $(BINDIR)/metalpdb2mol2.py

skip:
	@echo ""
	@echo "Skipping installation of pyMSMT."
	@echo ""

