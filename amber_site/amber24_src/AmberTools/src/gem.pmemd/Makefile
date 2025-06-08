SHELL=/bin/sh
CONFIG_FILE=../config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
CUDA_TOKEN=DCUDA

include $(CONFIG_FILE)

install: configured_serial $(NETCDF)
	$(MAKE) -C src/ install
	@echo "Installation of gem.pmemd complete"

parallel: configured_parallel $(NETCDF)
	$(MAKE) -C src/ parallel
	@echo "Installation of gem.pmemd.MPI complete"

uninstall.serial: clean
	-rm -f $(BINDIR)/gem.pmemd$(SFX)

uninstall.parallel: clean
	-rm -f $(BINDIR)/gem.pmemd.MPI$(SFX)

depends:
	$(MAKE) -C src/ depends

clean:
	$(MAKE) -C src/ clean

$(INCDIR)/netcdf.mod:
	$(MAKE) -C ../netcdf/src install

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
	    echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE)!" ;\
	    exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
	    { echo "Error: $(CONFIG_FILE) is not of type parallel!" ;\
	      echo "  Rerun $(CONFIG_COMMAND) and specify an MPI implementation." ;\
	      exit 2 ;\
	    } ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type parallel, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)
	@(if grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
	    echo "Error: $(CONFIG_FILE) is of type cuda, not serial!" ;\
	    echo "  Rerun $(CONFIG_COMMAND) and specify -nopar." ;\
	    exit 2 ;\
	fi ;\
	)
