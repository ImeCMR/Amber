#   top-level $(MAKE)file for the AmberTools suite of programs
AMBERTOOLS=AmberTools24
CONFIG_FILE=config.h
CONFIG_COMMAND=./configure
PARALLEL_TOKEN=DMPI
MIC_TOKEN=mmic
MIC2_TOKEN=DMIC2
MIC_OFFLOAD_TOKEN=DMIC_offload
PYINSTALL=PYTHONPATH=$(PYTHONPATH) $(PYTHON) setup.py install -f $(PYTHON_INSTALL)

include $(CONFIG_FILE)

.NOTPARALLEL:configured_serial
install: $(INSTALLTYPE)

sanderonly: configured_serial utils
	(cd pbsa && $(MAKE) -j1 libinstall )
	(cd cew && $(MAKE) install )
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)
	(cd rism && $(MAKE) librism )
	(cd sander && $(MAKE) $(AMBERHOME)/bin/sander)

sander.MPIonly: configured_parallel
	(cd pbsa && $(MAKE) -j1 libinstall )
	(cd cew && $(MAKE) install )
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)
	(cd rism && $(MAKE) librism_mpi )
	(cd sander && $(MAKE) $(AMBERHOME)/bin/sander.MPI)

utils: configured_serial
	@echo "Starting installation of ${AMBERTOOLS} serial at `date`".
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )

serial: configured_serial utils
	(cd pbsa && $(MAKE) -j1 install )
	(cd gbnsr6 && $(MAKE) install )
	(cd cifparse && $(MAKE) install )
	(cd cew && $(MAKE) install )

#	quick
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)

#   tcpb-cpp
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)

#	antechamber:
	(cd antechamber && $(MAKE) install )
	(cd sqm && $(MAKE) install )

#   miscellaneous:
	(cd reduce && $(MAKE) install )
	(if [ "$(HAS_BOOST)" = 'yes' -a "$(MOFT)" = 'yes' ]; then \
		cd moft && $(MAKE) install ;\
	fi ;\
	)

#   leap
	(cd leap && $(MAKE) install )

#	old Amber programs
#	(cd lib && $(MAKE) install)
	(cd addles && $(MAKE) install)
	(cd sander && $(MAKE) serial)
	(cd nmr_aux && $(MAKE) install)
	(cd nmode && $(MAKE) install)

#   cpptraj:
	$(MAKE) build_cpptraj
	(cd ambpdb && $(MAKE) install)

#   nab:
	(cd sff && $(MAKE) install )
	(cd pbsa && $(MAKE) -j1 libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd nab && $(MAKE) install )
#	(cd nss && $(MAKE) install )
	(cd etc && $(MAKE) install )

#   rism
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX) ;\
	fi ;\
	)

#   mdgx:
	(cd mdgx && $(MAKE) $(MDGX) )

#   xtalutil:
	(cd xtalutil && $(MAKE) install )

#   saxs and saxs_md:
	(cd saxs && $(MAKE) $(MAKE_SAXS) )

#	mm_pbsa
	(cd mm_pbsa && $(MAKE) install)

#	FEW
	(cd FEW && $(MAKE) install)

#	amberlite
	(cd amberlite && $(MAKE) install)

#	paramfit
	(cd paramfit && $(MAKE) install)

#	cphstats and cestats
	(cd cphstats && $(MAKE) install)

#	nfe-umbrella-slice
	(cd nfe-umbrella-slice && $(MAKE) install)

	@(if [ "$(PMEMD_GEM)" = 'yes' ]; then \
		cd gem.pmemd && $(MAKE) install ;\
	fi ;\
	)

#Python code
	(if [ "$(SKIP_PYTHON)" = 'no' ]; then \
		$(MAKE) python_serial ;\
	fi ;\
	)

	@echo ""
	@echo "Installation of ${AMBERTOOLS} serial is complete at `date`."
	@echo ""

# Serial install using cray compilers. Not working:
# pbsa, gbnsr6, sander, nmr_aux, nss, etc, mmpbsa, amberlite, quick
cray_serial: configured_serial
	@echo "Starting installation of ${AMBERTOOLS} serial (cray) at `date`".
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd cifparse && $(MAKE) install )
	(cd antechamber && $(MAKE) install )
	(cd sqm && $(MAKE) install )
	(cd cew && $(MAKE) install )
	(cd reduce && $(MAKE) install )
	(cd leap && $(MAKE) install )
	(cd lib && $(MAKE) install)
	(cd addles && $(MAKE) install)
	(cd nmode && $(MAKE) install)
	$(MAKE) build_cpptraj
	(cd ambpdb && $(MAKE) install)
	(cd sff && $(MAKE) install )
	(cd rism && $(MAKE) $(RISM) )
	(cd nab && $(MAKE) install )
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX) ;\
	fi ;\
	)
	(cd mdgx && $(MAKE) $(MDGX) )
	(cd xtalutil && $(MAKE) install )
	(cd saxs && $(MAKE) $(MAKE_SAXS) )
	(cd FEW && $(MAKE) install)
	(cd paramfit && $(MAKE) install)
	(cd cphstats && $(MAKE) install)
	(if [ "$(SKIP_PYTHON)" = 'no' ]; then \
		$(MAKE) python_serial ;\
	fi ;\
	)
	@echo ""
	@echo "Installation of ${AMBERTOOLS} serial (cray) is complete at `date`."
	@echo ""

cray_parallel: configured_parallel
	@echo "Starting installation of ${AMBERTOOLS} MPI (cray) at `date`".
	(cd etc && $(MAKE) parallel)
	(cd sff && $(MAKE) libsff_mpi )
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd cifparse && $(MAKE) install )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	$(MAKE) build_cpptraj
	(cd mdgx && $(MAKE) $(MDGX) )
	(cd nab && $(MAKE) install )
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) librism_mpi && $(MAKE) rism3d.snglpnt.MPI$(SFX) ;\
	fi ;\
	)
	(if [ "$(SKIP_PYTHON)" = 'no' ]; then \
		$(MAKE) python_parallel ;\
	fi ;\
	)
	@echo ""
	@echo "Installation of ${AMBERTOOLS} MPI (cray) is complete at `date`."
	@echo ""

cray_openmp: configured_serial
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd rism && $(MAKE) $(RISM) )
	(cd cifparse && $(MAKE) install )
	(cd sff && $(MAKE) install )
	(cd nab && $(MAKE) install )
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX) ;\
	fi ;\
	)
	$(MAKE) build_cpptraj
	$(MAKE) $(PYTRAJ)
	(cd saxs && $(MAKE) $(MAKE_SAXS) )
	(cd paramfit && $(MAKE) install )
	@echo ""
	@echo "Installation of ${AMBERTOOLS} OpenMP (cray) is complete at `date`."
	@echo ""

nabonly: configured_serial
#	utility routines and libraries:
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd c9x-complex && $(MAKE) $(C9XCOMPLEX) )
	(cd pbsa && $(MAKE) -j1 libinstall )
	(cd rism && $(MAKE) $(RISM) )
	(cd cifparse && $(MAKE) install )
	(cd sff && $(MAKE) install )
	(cd nab && $(MAKE) install )
#	(cd nss && $(MAKE) install )
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) rism3d.snglpnt$(SFX) ;\
	fi ;\
	)
	(cd leap && $(MAKE) install )
	(cd reduce && $(MAKE) install )

	@echo
	@echo "Completed installation of NAB"
	@echo

openmp: nabonly math_libraries
	$(MAKE) build_cpptraj
	-$(MAKE) $(PYTRAJ)
	(cd saxs && $(MAKE) $(MAKE_SAXS) )
	(cd paramfit && $(MAKE) install )

	@echo ""
	@echo "Installation of ${AMBERTOOLS} OpenMP is complete at `date`."
	@echo ""

.NOTPARALLEL:configured_parallel

parallel: configured_parallel mpinab
	@echo "Starting installation of ${AMBERTOOLS} MPI at `date`".
	(cd etc && $(MAKE) parallel)
	(cd pbsa && $(MAKE) -j1 libinstall )
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )
	(cd cew && $(MAKE) install )
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)
	(cd mdgx && $(MAKE) $(MDGX) )
	(cd addles && $(MAKE) install)
	(cd sander && $(MAKE) parallel)
	$(MAKE) build_cpptraj
	(if [ "$(RISM)" = 'yes' ]; then \
		cd rism && $(MAKE) librism_mpi && $(MAKE) rism3d.snglpnt.MPI$(SFX) ;\
	fi ;\
	)
	(if [ "$(SKIP_PYTHON)" = 'no' ]; then \
		$(MAKE) python_parallel ;\
	fi ;\
	)
	@(if [ "$(PMEMD_GEM)" = 'yes' ]; then \
		cd gem.pmemd && $(MAKE) parallel ;\
	fi ;\
	)
	@echo ""
	@echo "Installation of ${AMBERTOOLS} MPI is complete at `date`."
#	@echo "  (to install mpinab (if needed):"
#	@echo "   cd ${AMBERHOME}/AmberTools/src; make mpinab )"
	@echo ""

mpinab:
	(cd sff && $(MAKE) libsff_mpi )
	(cd ucpp-1.3 && $(MAKE) $(UCPP) )
	(cd cifparse && $(MAKE) install )
	(cd arpack && $(MAKE) install )
	(cd nab && $(MAKE) install )

configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\ # $(CONFIG_COMMAND) ;\
	fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) without the -mpi flag." ;\
		exit 2 ;\
	fi ;\
	)

configured_parallel: configured
	@(grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is of type serial, not parallel !" ;\
			echo "  Rerun $(CONFIG_COMMAND) with the -mpi flag." ;\
			exit 2 ;\
		} ;\
	)

configured_mic: configured
	@(grep $(MIC_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type mic_native !" ;\
			echo "Rerun $(CONFIG_FILE) WITH THE -mic_native flag." ;\
			exit 2 ;\
		} ;\
	)

configured_mic2: configured
	@(grep $(MIC2_TOKEN) $(CONFIG_FILE) > /dev/null || \
		{ echo "Error: $(CONFIG_FILE) is not of type mic2 !" ;\
			echo "Rerun $(CONFIG_FILE) WITH THE -mic2 flag." ;\
			exit 2 ;\
		} ;\
	)

configured_mic_offload: configured
	@(grep $(MIC_OFFLOAD_TOKEN) $(CONFIG_FILE) > /dev/null || \
                { echo "Error: $(CONFIG_FILE) is not of type mic_offload !" ;\
                        echo "Rerun $(CONFIG_FILE) WITH THE -mic_offload flag." ;\
                        exit 2 ;\
                } ;\
        )

netcdf_uninstall:
	-(cd netcdf-4.6.1 && $(MAKE) uninstall)
	-(cd netcdf-fortran-4.4.4 && $(MAKE) uninstall)
	-/bin/rm -f netcdf_config.log ../lib/libnetcdf* ../include/netcdf.* \
	            $(BINDIR)/nf-config

clean::
	-(cd ucpp-1.3 && $(MAKE) clean )
	-(cd cifparse && $(MAKE) clean )
	-(cd nab && $(MAKE) clean )
	-(cd sff && $(MAKE) clean )
#	-(cd nss && $(MAKE) clean )
	-(cd leap && $(MAKE) clean )
	-(cd arpack && $(MAKE) clean )
	-(cd blas && $(MAKE) clean )
	-(cd xblas && $(MAKE) clean )
	-(cd lapack && $(MAKE) clean )
	-(cd c9x-complex && $(MAKE) clean )
	-(cd etc && $(MAKE) clean )
	-(cd chamber && $(MAKE) clean )
	-(cd pbsa && $(MAKE) clean )
	-(cd gbnsr6 && $(MAKE) clean )
	-(cd antechamber && $(MAKE) clean )
	-(cd sqm && $(MAKE) clean )
	-(cd cew && $(MAKE) clean )
	-(cd lib && $(MAKE) clean )
	-(cd ptraj && $(MAKE) clean )
	-(cd cpptraj && $(MAKE) clean)
	-(cd ambpdb && $(MAKE) clean)
	-(cd reduce && $(MAKE) clean )
	-(cd mdgx && $(MAKE) clean )
	-(cd xtalutil && $(MAKE) clean )
	-(cd saxs && $(MAKE) clean )
	-(cd amberlite && /bin/rm -fr build/)
	-(cd paramfit && $(MAKE) clean )
	-(cd rism && $(MAKE) clean )
	-(cd fftw-3.3 && $(MAKE) clean )
	-($(RM) -f fftw-3.3/mpi/fftw3-mpi.f03 )
	-(cd mmpbsa_py && /bin/rm -fr MMPBSA.py.MPI build/ mpi4py-1.2.2/)
	-(cd mm_pbsa && $(MAKE) clean )
	-(cd FEW && $(MAKE) clean )
	-(cd moft && $(MAKE) clean )
	-(cd packmol_memgen && /bin/rm -fr build/ build.log && cd packmol_memgen && make clean)
	-(cd parmed && /bin/rm -fr build/ build.log && $(PYTHON) setup.py clean)
	-(cd pysander && /bin/rm -fr build/ build.log )
	-(/bin/rm -rf pytraj/scripts/*.pyc)
	-(cd pytraj && /bin/rm -fr build/ build.log)
	-(cd pdb4amber && /bin/rm -fr build/ build.log)
	-(cd cphstats && $(MAKE) clean )
	-(cd lib && $(MAKE) clean)
	-(cd addles && $(MAKE) clean)
	-(cd sander && $(MAKE) clean)
	-(cd nmr_aux && $(MAKE) clean)
	-(cd nmode && $(MAKE) clean)
	-(cd emil && $(MAKE) clean)
	-(cd sebomd && $(MAKE) clean)
	-(cd cew && $(MAKE) clean )
	-(if [ "$(QUICK)" = 'yes' ]; then \
			cd quick && $(MAKE) clean ;\
		fi ;\
		)
	-(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) clean ;\
	    fi ;\
	    )
	-(cd nfe-umbrella-slice && $(MAKE) clean)
	-(cd pymsmt && $(MAKE) clean)
	-(cd netcdf-4.6.1 && $(MAKE) clean)
	-(cd netcdf-fortran-4.4.4 && $(MAKE) clean)
	-(cd gem.pmemd && $(MAKE) clean)

silentclean:
	@$(MAKE) clean 1> /dev/null 2>&1

#  uninstall Amber stuff, but leave behind (for example) files
#     created by configure_mpich, configure_python, configure_openmpi
uninstall: netcdf_uninstall
	-(cd ucpp-1.3 && $(MAKE) uninstall )
	-(cd cifparse && $(MAKE) uninstall )
	-(cd nab && $(MAKE) uninstall )
	-(cd sff && $(MAKE) uninstall )
#	-(cd nss && $(MAKE) uninstall )
	-(cd leap && $(MAKE) uninstall )
	-(cd arpack && $(MAKE) uninstall )
	-(cd blas && $(MAKE) uninstall )
	-($(RM) -f $(LIBDIR)/libxblas.a )
	-(cd lapack && $(MAKE) uninstall )
	-(cd c9x-complex && $(MAKE) uninstall )
	-(cd etc && $(MAKE) uninstall )
	-(cd chamber && $(MAKE) uninstall )
	-(cd pbsa && $(MAKE) uninstall )
	-(cd gbnsr6 && $(MAKE) uninstall )
	-(cd antechamber && $(MAKE) uninstall )
	-(cd sqm && $(MAKE) uninstall )
	-(cd cew && $(MAKE) uninstall )
	-(cd lib && $(MAKE) uninstall )
	-(cd ptraj && $(MAKE) uninstall )
	-(cd reduce && $(MAKE) uninstall )
	-(cd mdgx && $(MAKE) uninstall )
	-(cd xtalutil && $(MAKE) uninstall )
	-(cd saxs && $(MAKE) uninstall )
	-(cd paramfit && $(MAKE) uninstall )
	-(cd rism && $(MAKE) uninstall )
	-(cd fftw-3.3 && $(MAKE) uninstall)
	-(cd cpptraj && $(MAKE) uninstall)
	-(cd ambpdb && $(MAKE) uninstall)
	-(cd mmpbsa_py && $(MAKE) uninstall)
	-(cd mm_pbsa && $(MAKE) uninstall)
	-(cd FEW && $(MAKE) uninstall)
	-(cd moft && $(MAKE) uninstall)
	-(/bin/rm -fr $(BINDIR)/*parmed $(BINDIR)/pdb4amber $(BINDIR)/packmol-memgen \
		$(BINDIR)/packmol $(BINDIR)/memembed)
	-(cd amberlite && $(MAKE) uninstall )
	-/bin/rm -f ../test/numprocs
	-/bin/rm -fr $(BINDIR)/ncmpidump $(BINDIR)/ncmpigen $(BINDIR)/ncvalid $(LIBDIR)/libnetcdf*
	-(cd $(BINDIR) && /bin/rm -f charmmlipid2amber.py cpinutil.py ceinutil.py \
		cpeinutil.py mdout_analyzer.py finddgref.py fitpkaeo.py \
		fixremdcouts.py genremdinputs.py softcore_setup.py )
	-(cd $(LIBDIR) && /bin/rm -f libxblas-amb.a libfftw3_mpi.a libfftw3_mpi.la)
	-(cd $(LIBDIR) && /bin/rm -f libmdgx.so)
	-(cd $(LIBDIR) && /bin/rm -f libemil.a)
	-(/bin/rm -f $(LIBDIR)/libcpptraj$(SHARED_SUFFIX))
	-(cd $(INCDIR) && /bin/rm -f fftw3-mpi.f03 fftw3-mpi.h)
	-(cd lib && $(MAKE) uninstall)
	-(cd addles && $(MAKE) uninstall)
	-(cd sander && $(MAKE) uninstall.serial uninstall.cuda uninstall.parallel)
	-(cd cew && $(MAKE) uninstall )
	-(if [ "$(QUICK)" = 'yes' ]; then \
			cd quick && $(MAKE) uninstall ;\
		fi ;\
		)
	-(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) uninstall ;\
	    fi ;\
	    )
	-(cd nfe-umbrella-slice && $(MAKE) uninstall)
	-(cd nmr_aux && $(MAKE) uninstall)
	-(cd nmode && $(MAKE) uninstall)
	-(cd cphstats && $(MAKE) uninstall)
	-(cd pymsmt && $(MAKE) uninstall)
	-(cd gem.pmemd && $(MAKE) uninstall)
	-(cd $(LIBDIR) && /bin/rm -fr python*/)
	-/bin/rm -f config.h ../../config.h

math_libraries:
	(cd arpack && $(MAKE) install )
	(cd lapack && $(MAKE) $(LAPACK) )
	(cd blas && $(MAKE) $(BLAS) )

cuda_serial: math_libraries
	(cd pbsa && $(MAKE) clean && $(MAKE) -j1 install )
	(cd mdgx && $(MAKE) clean && $(MAKE) -j1 cuda )
	$(MAKE) build_cpptraj
	(cd cew && $(MAKE) install )
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd sander && $(MAKE) cuda_serial ;\
	fi ;\
	)
	@echo "AmberTools serial CUDA install complete."

cuda_parallel: math_libraries
	$(MAKE) build_cpptraj
	(cd cew && $(MAKE) install )
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd quick && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(TCPB)" = 'yes' ]; then \
		cd tcpb-cpp && $(MAKE) install ;\
	fi ;\
	)
	(if [ "$(QUICK)" = 'yes' ]; then \
		cd sander && $(MAKE) cuda_parallel ;\
	fi ;\
	)
	@echo "AmberTools parallel CUDA install complete."

ompmpi: math_libraries
	$(MAKE) build_cpptraj
	@echo "AmberTools OpenMP/MPI install complete."

mic2: configured_mic2
	@echo "$(AMBERTOOLS) has no mic2-enabled components"

mic: configured_mic
	@echo "$(AMBERTOOLS) has no mic-enabled components"

mic_parallel: configured_mic configured_parallel
	@echo "$(AMBERTOOLS) has no mic-enabled components"

mic_offload: configured_mic_offload
	@echo "$(AMBERTOOLS) has no mic-enabled components"

#Rule for building python code
#all_python: parmed pymsmt pytraj pysander

python_serial:: parmed mmpbsa_py pysander pymsmt $(PYTRAJ) pdb4amber packmol-memgen
	$(PYTHON) fix_conda_fortran.py

python_parallel:: parmed mmpbsa_py_parallel

# Rule for building cpptraj.
build_cpptraj:
	(cd cpptraj && $(MAKE) install)

# Rule for building libcpptraj
libcpptraj:
	(cd cpptraj && $(MAKE) libcpptraj)

# Add a rule to build chamber for anyone that wants to compile it easily
chamber::
	(cd chamber && $(MAKE) install )

# Python package standalone rules
pytraj::
	$(MAKE) $(LIBCPPTRAJ)
ifeq ($(INSTALLTYPE), openmp)
	(cd pytraj && AMBERHOME=$(AMBER_SOURCE) AMBER_PREFIX=$(AMBER_PREFIX) $(PYINSTALL) $(PYTHONLOG)
else ifeq ($(INSTALLTYPE), serial)
	(cd pytraj && AMBERHOME=$(AMBER_SOURCE) AMBER_PREFIX=$(AMBER_PREFIX) $(PYINSTALL) --disable-openmp $(PYTHONLOG) )
endif

no_pytraj:
	@echo "No python development headers installed; skipping pytraj."

no_libcpptraj:
	@echo "Pytraj is not being installed. Skipping libcpptraj."
packmol-memgen::
	(if [ "$(HAS_BOOST)" = 'yes' ]; then \
		cd packmol_memgen && $(PYINSTALL) $(PYTHONLOG) && cp packmol_memgen/lib/packmol/packmol packmol_memgen/lib/memembed/bin/memembed $(BINDIR) ;\
	else \
		cd packmol_memgen && $(PYINSTALL) --no-memembed $(PYTHONLOG) && cp packmol_memgen/lib/packmol/packmol $(BINDIR) ;\
	fi ;\
	)
parmed::
	(cd parmed && $(PYINSTALL) $(PYTHONLOG))
pymsmt::
	(cd pymsmt && $(PYINSTALL) $(PYTHONLOG) )
pysander::
	(cd pysander && $(MAKE) $(PYSANDER) $(PYTHONLOG))
mmpbsa_py::
	(cd mmpbsa_py && $(MAKE) install)
mmpbsa_py_parallel::
	(cd mmpbsa_py && $(MAKE) parallel)
pdb4amber::
	(cd pdb4amber && $(PYINSTALL) $(PYTHONLOG))
pdb4amber.develop::
	# Install development version so you can edit file in source folder, run code without install again.
	# Note: if you install develop version, do not use "pip uninstall", you must use
	# make pdb4amber.develop_uninstall
	(cd pdb4amber && PYTHONPATH=$(PYTHONPATH) $(PYTHON) setup.py develop $(PYTHON_INSTALL))
pdb4amber.develop_uninstall:
	(cd pdb4amber && PYTHONPATH=$(PYTHONPATH) $(PYTHON) setup.py develop $(PYTHON_INSTALL) --uninstall)
