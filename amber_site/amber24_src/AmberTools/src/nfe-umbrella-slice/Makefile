# Makefile for nfe-umbrella-slice
include ../config.h

target = nfe-umbrella-slice$(SFX)
nfe_umbrella_slice_OBJECTS = umbrella.o \
	dimension.o utils.o main.o

DEFS = -DHAVE_CONFIG_H
INCLUDES = -I. $(NETCDFINC)
CXXLINK = $(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@
CXXCOMPILE = $(CXX) $(DEFS) $(INCLUDES) $(CXXFLAGS)

.SUFFIXES: .cc .o
.PRECIOUS: Makefile

all: install

install: $(target)

nfe-umbrella-slice$(SFX): $(nfe_umbrella_slice_OBJECTS) 
	@echo "[NFE]  CXX $@"
	$(VB)$(CXXLINK) $(nfe_umbrella_slice_OBJECTS) $(NETCDFLIB)
	$(VB)$(MV) nfe-umbrella-slice$(SFX) $(BINDIR)

.cc.o:
	@echo "[NFE]  CXX $<"
	$(VB)$(CXXCOMPILE) -c -o $@ $<
	
clean:
	$(RM) -f *.o
	
uninstall: clean
	$(RM) -f $(BINDIR)/nfe-umbrella-slice$(SFX)
	
.PHONY: all install clean uninstall	
