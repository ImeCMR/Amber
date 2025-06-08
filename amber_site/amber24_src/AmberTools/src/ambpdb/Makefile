# AmbPDB Makefile
include ../config.h

AMBPDBSOURCES=AmbPDB.cpp Trajout_Single_Stdout.cpp DSL.cpp

AMBPDB_OBJECTS=$(AMBPDBSOURCES:.cpp=.o)

CPPTRAJ_BASE=$(BASEDIR)/AmberTools/src/cpptraj

CPPTRAJ_INCDIR=$(CPPTRAJ_BASE)/src

# Static cpptraj libraries ambpdb depends on
LIBCPPTRAJ= \
  $(LIBDIR)/libcpptraj_parm.a \
  $(LIBDIR)/libcpptraj_traj.a \
  $(LIBDIR)/libcpptraj_file.a \
  $(LIBDIR)/libcpptraj_core.a


# Default target
all: install

install: install_ambpdb

# Cpptraj static libraries -------------
libcpptraj.static:
	cd $(CPPTRAJ_BASE) && $(MAKE) libstatic

# AmbPDB -------------------------------
install_ambpdb: ambpdb
	mv ambpdb$(SFX) $(BINDIR)

ambpdb: $(AMBPDB_OBJECTS) libcpptraj.static
	@echo "[AMBPDB]  CXX $@"
	$(VB)$(CXX) -o ambpdb$(SFX) $(AMBPDB_OBJECTS) $(LIBCPPTRAJ) $(NETCDFLIB) $(ZLIB) $(BZLIB) $(LDFLAGS)

# Dependencies
AmbPDB.o: AmbPDB.cpp Trajout_Single_Stdout.h

Trajout_Single_Stdout.o: Trajout_Single_Stdout.cpp Trajout_Single_Stdout.h DSL.h

DSL.o: DSL.cpp DSL.h

# Build rules
#.cpp.o:
#	@echo CXX $<
#	@$(CXX) $(CXXFLAGS) -c -o $@ $<
.cpp.o:
	@echo "[AMBPDB]  CXX $<"
	$(VB)$(CXX) -I$(CPPTRAJ_INCDIR) $(CXXOPTFLAGS) $(CXXFLAGS) -c -o $@ $<

# Clean/uninstall targets
clean:
	/bin/rm -f $(AMBPDB_OBJECTS) ambpdb$(SFX) 

uninstall: clean
	/bin/rm -f $(BINDIR)/ambpdb$(SFX)

