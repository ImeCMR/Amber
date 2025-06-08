include ../config.h

install: all
	@/bin/mv saxs_rism${SFX} saxs_md${SFX} $(BINDIR)
	@/bin/mv libsaxs_rism.so $(LIBDIR)
	@/bin/mv libsaxs_md.so $(LIBDIR)

no-shared: saxs_rism$(SFX) saxs_md$(SFX)
	/bin/mv saxs_rism${SFX} saxs_md${SFX} $(BINDIR)

lib: libsaxs_md.so
	@/bin/mv libsaxs_md.so $(LIBDIR)

all: saxs_rism$(SFX) saxs_md$(SFX) libsaxs_rism.so libsaxs_md.so

saxs_rism$(SFX): saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
	   -o $@ saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o

saxs_md$(SFX): saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
	   -o $@ saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o

libsaxs_md.so: saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) -shared $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
		-o $@ saxs_md.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o

libsaxs_rism.so: saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o
	@echo "[SAXS]  CXX $@"
	$(VB)$(CXX) -shared $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) $(LDFLAGS) \
		-o $@ saxs_rism.o lebedev.o atm_f.o dx.o f_calc.o IO.o pdb.o vol_decomp.o

.cpp.o:
	@echo "[SAXS]  CXX $<"
	$(VB)$(CXX) -c $(CXXOPTFLAGS) $(CXXFLAGS) $(AMBERCXXFLAGS) -I$(INCDIR) -o $@ $<
	
#lebedev.cpp:  lebedev.hpp

pysaxs: lib
	(cd interface && python ./setup.py install --prefix=$(AMBERHOME))

clean:
	/bin/rm -f *.o
	/bin/rm -rf interface/build

uninstall:
	/bin/rm -f $(BINDIR)/saxs_rism$(SFX) $(BINDIR)/saxs_md$(SFX)
	/bin/rm -f $(LIBDIR)/libsaxs_md.so
	/bin/rm -f $(LIBDIR)/libsaxs_rism.so
