#
#************************************************************************
#                              AMBER                                   **
#                                                                      **
#               Copyright (c) 1986, 1991, 1995, 1997, 1999             **
#                Regents of the University of California               **
#                       All Rights Reserved.                           **
#                                                                      **
#  This software provided pursuant to a license agreement containing   **
#  restrictions on its disclosure, duplication, and use. This software **
#  contains confidential and proprietary information, and may not be   **
#  extracted or distributed, in whole or in part, for any purpose      **
#  whatsoever, without the express written permission of the authors.  **
#  This notice, and the associated author list, must be attached to    **
#  all copies, or extracts, of this software. Any additional           **
#  restrictions set forth in the license agreement also apply to this  **
#  software.                                                           **
#************************************************************************
#
# Makefile for SQM (stand-alone quantum program)
#
include ../config.h

# Sqm Fortran source files are free format
LOCALFLAGS = 

SHELL=/bin/sh
VPATH=

CPPFLAGS += -I. -I$(INCDIR)


#LIBCEWOBJ = Objects specific to building as a library.
LIBCEWOBJ = Ewald$(OSFX) iface$(OSFX)



ifneq (, $(filter $(INSTALLTYPE), parallel cuda_parallel))
LSFX=_mpi
CPPFLAGS += -DWITH_MPI
OBJS=$(LIBCEWOBJ)
else
LSFX=
OBJS=$(LIBCEWOBJ)
endif



install: $(LIBDIR)/libcew$(LSFX)$(SHARED_SUFFIX)


$(LIBDIR)/libcew$(LSFX)$(SHARED_SUFFIX): $(OBJS)
	@echo "[CEW]  AR $@"
	$(CXX) $(MAKE_SHARED) $(filter-out -DMPI, $(CXXOPTFLAGS) $(CXXFLAGS)) \
	-o $@ $^ \
	-L$(LIBDIR) $(FLIBS_FFTW3) $(FLIBSF) $(NOFORTRANMAIN) -lrt


%$(OSFX): CXXFLAGS += $(CXXOPTFLAGS) 


%_mpi.o: %.cpp
	$(CXX) $(filter-out -DMPI, $(CPPFLAGS) $(CXXFLAGS)) -o $@ -c $<


clean:
	/bin/rm -f *.$(OSFX)

uninstall:
	-rm -f $(LIBDIR)/libcew$(LSFX)$(SHARED_SUFFIX)



