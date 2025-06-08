#
#************************************************************************
#
#************************************************************************
#
include ../config.h

CONFIG_FILE=../config.h
CONFIG_COMMAND=configure
PARALLEL_TOKEN=DMPI

CUDA_TOKEN=DCUDA
MATRIX_FMT=DDIA#DDIA/DCSR

ifneq (,$(findstring hipcc,$(NVCC)))
  CUSP_LIB = -I../cusplibrary-hip # select this if hipcc
  CUDA_INC=-I$(HIP_HOME)/include
else
  CUSP_LIB = -I../cusplibrary-cuda9
  CUDA_INC=-I$(CUDA_HOME)/include
endif

# ---- Begin file definitions ----
# PBSA Fortran source files are free format
LOCALFLAGS = $(FREEFORMAT_FLAG)

SHELL=/bin/sh
VPATH=

.SUFFIXES:
.SUFFIXES: .F90 .o

OBJ=\
	dsvdc.o  pbsa.o pbsa_lib.o pb_read.o pb_write.o getcoor.o \
	runmd.o runmin.o force.o pb_init.o sa_driver.o density.o np_force.o pb_force.o \
        pb_p3m.o pb_fdfrc.o pb_crgview.o \
        pb_fddrv.o pb_anaiim.o pb_nhadrv.o pb_iimdrv.o pb_augdrv.o pb_nlsolver.o \
	pb_bldsys.o pb_mpfrc.o pb_direct.o pb_list.o timer.o \
	egb.o ene.o locmem.o myopen.o rdparm.o pb_lsolver.o \
	debug.o \
	irre31.o irre32.o  phidev.o  project.o  qld.o miniop.o \
	transf.o curv.o GrToPr.o jumps.o \
	regular.o bicg.o gmres.o \
	daxpy.o ddot.o dheqr.o dpigmr.o dxlcal.o \
	dchkw.o dgmres.o dnrm2.o drlcal.o isdgmr.o \
	dcopy.o dhels.o dorth.o \
        ds2y.o dscal.o dsilus.o dslugm.o dslui.o dslui2.o dsmv.o qs2i1d.o \
	d1mach.o fdump.o i1mach.o j4save.o xercnt.o xerhlt.o xermsg.o  \
        xerprn.o xersve.o xgetua.o \
        pb_chkfrc.o svdcmp.o svbksb.o pythag.o\
        interpX.o gen_dx_file.o \
        parms.o memory_module.o rdpqr.o \
	mlses_wrapper.o mlses.o utils.o

CUDAOBJ=\
	cuda_pb.o kLinearSolvers.o cuda_LinearSolvers.o cusparse_LinearSolvers.o cusp_LinearSolvers.o

# GM: when updating this variable, please don't omit to update makedepend
LIBPBSAOBJS= \
        interpX.LIBPBSA.o parms.LIBPBSA.o pb_fdfrc.LIBPBSA.o pb_force.LIBPBSA.o pb_p3m.LIBPBSA.o \
        pb_read.LIBPBSA.o sa_driver.LIBPBSA.o \
        curv.LIBPBSA.o density.LIBPBSA.o dsvdc.LIBPBSA.o gen_dx_file.LIBPBSA.o \
        GrToPr.LIBPBSA.o interface.o memory_module.LIBPBSA.o miniop.LIBPBSA.o myopen.LIBPBSA.o \
        np_force.LIBPBSA.o pb_bldsys.LIBPBSA.o pb_crgview.LIBPBSA.o pb_direct.LIBPBSA.o pb_fddrv.LIBPBSA.o \
        pb_init.LIBPBSA.o pb_list.LIBPBSA.o pb_lsolver.LIBPBSA.o pb_mpfrc.LIBPBSA.o pb_nlsolver.LIBPBSA.o  \
        pbsaapi.LIBPBSA.o pbsa_lib.LIBPBSA.o pb_write.LIBPBSA.o phidev.LIBPBSA.o project.LIBPBSA.o \
        pythag.LIBPBSA.o svbksb.LIBPBSA.o svdcmp.LIBPBSA.o \
        timer.LIBPBSA.o transf.LIBPBSA.o

# GM: when updating this variable, please don't omit to update makedepend
SANDERPBSAOBJS= \
        interpX.SANDER.o np_force.SANDER.o  pb_fdfrc.SANDER.o pb_force.SANDER.o pb_init.SANDER.o \
        pb_p3m.SANDER.o pb_read.SANDER.o pb_write.SANDER.o sa_driver.SANDER.o timer.SANDER.o \
        curv.SANDER.o density.SANDER.o dsvdc.SANDER.o gen_dx_file.SANDER.o GrToPr.SANDER.o \
        memory_module.SANDER.o miniop.SANDER.o pb_bldsys.SANDER.o pb_crgview.SANDER.o \
        pb_direct.SANDER.o pb_fddrv.SANDER.o  pb_list.SANDER.o pb_lsolver.SANDER.o  pb_mpfrc.SANDER.o \
        pb_nlsolver.SANDER.o pbsa_lib.SANDER.o phidev.SANDER.o project.SANDER.o pythag.SANDER.o  \
        svbksb.SANDER.o svdcmp.SANDER.o transf.SANDER.o
# ---- End file definitions ----

# ---- Begin installation options ----
install: install.$(INSTALLTYPE)

install.serial: serial
	/bin/mv pbsa$(SFX) $(BINDIR)

install.cuda_serial: cuda
	/bin/mv pbsa.cuda$(SFX) $(BINDIR)

serial: configured_serial pbsa$(SFX)

cuda: configured_cuda pbsa.cuda$(SFX)

uninstall:
	-/bin/rm -f $(BINDIR)/pbsa$(SFX) pbsa.cuda$(SFX) $(LIBDIR)/libpbsa.a $(LIBDIR)/libFpbsa.a

libinstall: libinstall.$(INSTALLTYPE)

libinstall.openmp: libinstall.serial

libinstall.serial: libpbsa.a libFpbsa.a
	/bin/mv libpbsa.a libFpbsa.a $(LIBDIR)

libinstall.parallel: libpbsa.a libFpbsa.a
	/bin/mv libpbsa.a libFpbsa.a $(LIBDIR)
# ---- End installation options ----

# ---- Begin configuration check ----
configured:
	@(if [ ! -f $(CONFIG_FILE) ] ; then \
		echo "Error: $(CONFIG_COMMAND) must be executed before $(MAKE) !" ;\
		exit 2 ;\
		fi ;\
	)

configured_serial: configured
	@(if grep $(PARALLEL_TOKEN) $(CONFIG_FILE) > /dev/null ; then \
		echo "Error: $(CONFIG_FILE) is of type parallel, not serial !" ;\
		echo "  Rerun $(CONFIG_COMMAND) and do NOT specify -mpi." ;\
		exit 2 ;\
	fi ;\
	)

configured_cuda: configured_serial
ifneq (,$(findstring hipcc,$(NVCC)))
			@(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
	                { echo "Error: $(CONFIG_FILE) is not of type hip !" ;\
	                        echo " Rerun $(CONFIG_COMMAND) and specify -hip." ;\
	                        exit 2 ;\
	                };\
	        if [ -z $(HIP_HOME) ]; then \
	                echo "Error: HIP_HOME is not set. This must point to your AMD ROCm installation"; \
	                exit 2; \
	                fi; \
	        )
else
	        @(grep $(CUDA_TOKEN) $(CONFIG_FILE) > /dev/null || \
	                { echo "Error: $(CONFIG_FILE) is not of type cuda !" ;\
	                        echo " Rerun $(CONFIG_COMMAND) and specify -cuda." ;\
	                        exit 2 ;\
	                };\
	        if [ -z $(CUDA_HOME) ]; then \
	                echo "Error: CUDA_HOME is not set. This must point to your NVIDIA tools installation"; \
	                exit 2; \
	                fi; \
	        )
endif
# ---- End configuration check ----

# ---- Begin cuda object files ----
.UPDATE:
	$(info Always rebuild [ pb_lsolver.o pb_fddrv.o [ $(CUDAOBJ) ]] to make predefined macros effective)

$(CUDAOBJ) pb_fddrv.o pb_lsolver.o: .UPDATE

ifneq (,$(findstring hipcc,$(NVCC)))
  NVCCEXTRAOPTS = -I/opt/rocm/include -fPIC -DAMBER_PLATFORM_AMD
else
  NVCCEXTRAOPTS = -ccbin $(CXX)
endif

cuda_pb.o: cuda_pb.cu
	$(NVCC)  $(NVCCEXTRAOPTS) -o $@ -c $<
kLinearSolvers.o: kLinearSolvers.cu
	$(NVCC)  $(NVCCEXTRAOPTS) -o $@ -c $<
cusparse_LinearSolvers.o: cusparse_LinearSolvers.cu
	$(NVCC)  $(NVCCEXTRAOPTS) -o $@ -c $<
cusp_LinearSolvers.o: cusp_LinearSolvers.cu
	$(NVCC)  $(NVCCEXTRAOPTS) $(CUSP_LIB) -o $@ -c $< -$(MATRIX_FMT)
cuda_LinearSolvers.o: cuda_LinearSolvers.cu
	$(NVCC)  $(NVCCEXTRAOPTS) $(CUSP_LIB) -o $@ -c $<
# ---- End cuda object files ----

# ---- Begin build pbsa options ----
pbsa$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial libtorch
	@echo "[PBSA]  FC $@"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa$(SFX) $(OBJ) \
	        ../lib/nxtsec.o ../lib/random.o \
	        -L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)

pbsa.cuda$(SFX): $(OBJ) syslib netlib c9x-complex configured_serial $(CUDAOBJ)
	@echo "[PBSA]  FC $@"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o pbsa.cuda$(SFX) \
		$(OBJ) $(CUDAOBJ) $(PBSA_CU_LIBS) $(CUDA_INC) \
		../lib/nxtsec.o ../lib/random.o \
		-L$(LIBDIR) $(FLIBSF) $(LDFLAGS) $(AMBERLDFLAGS)
# ---- End build pbsa options ----

# ---- Begin build other options ----
simplepbsa$(SFX): simplepbsa.o gopt.o libpbsa.a sfflib syslib netlib c9x-complex
	@echo "[PBSA]  FC $@"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o simplepbsa$(SFX) simplepbsa.o gopt.o \
		libpbsa.a $(LIBDIR)/libsff.a ../lib/nxtsec.o -L$(LIBDIR) $(FLIBSF) \
		$(LDFLAGS) $(AMBERLDFLAGS)
	/bin/mv simplepbsa$(SFX) $(BINDIR)

libpbsa.a: $(LIBPBSAOBJS) decomp.o
	@echo "[PBSA]  AR $@"
	$(VB)rm -f libpbsa.a
	$(VB)$(AR) libpbsa.a $(LIBPBSAOBJS) decomp.o > /dev/null
	$(VB)$(RANLIB) libpbsa.a

#  Following should only be invoked from the sander Makefile:
libFpbsa.a: $(SANDERPBSAOBJS) decomp.o file_io_dat.o memory_module.SANDER.o
	@echo "[PBSA]  AR $@"
	$(VB)rm -f libFpbsa.a
	$(VB)$(AR) libFpbsa.a $(SANDERPBSAOBJS) decomp.o file_io_dat.o memory_module.SANDER.o > /dev/null
	$(VB)$(RANLIB) libFpbsa.a

# ---- End build other options ----

syslib:
	cd ../lib && $(MAKE) nxtsec.o random.o

sfflib:
	cd ../sff && $(MAKE) install

netlib:
	cd ../lapack && $(MAKE) $(LAPACK)
	cd ../blas && $(MAKE) $(BLAS)
	cd ../arpack && $(MAKE) install

c9x-complex:
	@if test $(C9XCOMPLEX) != "skip"; then \
	   cd ../c9x-complex && $(MAKE) libmc.a; \
	 fi

# --- Predefine -DCUDA for pb_lsolver.F90/pb_fddrv.F90 ---
ifeq ($(INSTALLTYPE), cuda_serial)
    CU_define=-$(CUDA_TOKEN)
else
    CU_define=
endif

pb_init.o: pb_init.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pb_init.SANDER.o: pb_init.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DSANDER -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

pb_init.LIBPBSA.o: pb_init.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

np_force.SANDER.o: np_force.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DSANDER -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

np_force.LIBPBSA.o: np_force.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c -O0 $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

%.LIBPBSA.o: %.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DLIBPBSA -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $<

%.SANDER.o: %.F90
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -DSANDER -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -I../sander \
	-o $@ $<

interface.o: interface.c
	@echo "[PBSA]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

mlses_wrapper.o: mlses_wrapper.cpp
	@echo "[PBSA]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

mlses.o: mlses.cpp
	@echo "[PBSA]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

utils.o: utils.cpp
	@echo "[PBSA]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

.cpp.o:   $<
	@echo "[PBSA]  CC $<"
	$(VB)$(CC) -c $(COPTFLAGS) $(CFLAGS) $(AMBERFFLAGS) -o $@ $<

.F90.o:   $<
	@echo "[PBSA]  FC $<"
	$(VB)$(FC) $(PBSAFLAG) $(FPPFLAGS) -c $(FOPTFLAGS) $(FFLAGS) $(AMBERFFLAGS) -o $@ $< $(CU_define) -$(MATRIX_FMT)

depend::
	./makedepend > depend

clean:
	/bin/rm -f pbsa$(SFX) libpbsa.a
	/bin/rm -f *.o
	/bin/rm -f *.mod *__genmod.F90
	/bin/rm -f *.d

distclean: clean



#----------HEADER DEPENDENCIES (currently from makedepend)
include depend
# DO NOT DELETE
