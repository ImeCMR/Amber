project(pbsa C CXX Fortran)

set(EXECUTABLE_SOURCES dsvdc.F90 pbsa.F90 pb_read.F90 pb_write.F90 getcoor.F90
	runmd.F90 runmin.F90 force.F90 pbsa_lib.F90 pb_init.F90 density.F90 np_force.F90 pb_force.F90
	pb_p3m.F90 pb_fdfrc.F90 pb_crgview.F90 pb_fddrv.F90
        pb_anaiim.F90 pb_nhadrv.F90 pb_iimdrv.F90 pb_augdrv.F90 pb_nlsolver.F90
	pb_bldsys.F90 pb_mpfrc.F90 pb_direct.F90 pb_list.F90 timer.F90
	egb.F90 ene.F90 locmem.F90 myopen.F90 rdparm.F90 pb_lsolver.F90
	debug.F90 irre32.F90  phidev.F90  project.F90  qld.F90
	irre31.F90  miniop.F90 transf.F90 curv.F90 GrToPr.F90 jumps.F90
	regular.F90 bicg.F90 gmres.F90 sa_driver.F90
	daxpy.F90  ddot.F90    dheqr.F90  dpigmr.F90  dxlcal.F90
	dchkw.F90  dgmres.F90  dnrm2.F90  drlcal.F90  isdgmr.F90
	dcopy.F90 dhels.F90 dorth.F90
        ds2y.F90 dscal.F90 dsilus.F90 dslugm.F90 dslui.F90 dslui2.F90 dsmv.F90 qs2i1d.F90
	d1mach.F90 fdump.F90 i1mach.F90 j4save.F90 xercnt.F90 xerhlt.F90 xermsg.F90
        xerprn.F90 xersve.F90 xgetua.F90
	pb_chkfrc.F90 svdcmp.F90 svbksb.F90 pythag.F90
	interpX.F90 gen_dx_file.F90
	parms.F90 memory_module.F90 rdpqr.F90)

#sources for libpbsa
set(LIBPBSA_SOURCES interpX.F90 parms.F90 pb_fdfrc.F90 pb_force.F90 pb_p3m.F90
	pb_read.F90 sa_driver.F90
	curv.F90 density.F90 dsvdc.F90 gen_dx_file.F90
	GrToPr.F90 interface.c memory_module.F90 miniop.F90 myopen.F90
	np_force.F90 pb_bldsys.F90 pb_crgview.F90 pb_direct.F90 pb_fddrv.F90
	pb_init.F90 pb_list.F90 pb_lsolver.F90 pb_mpfrc.F90 pb_nlsolver.F90
	pbsaapi.F90 pbsa_lib.F90 pb_write.F90 phidev.F90 project.F90
	pythag.F90 svbksb.F90 svdcmp.F90
	timer.F90 transf.F90)

#sources for the sander PBSA library (libFpbsa in the makefile system).
#will be built with SANDER defined.
set(SANDER_PBSA_SOURCES interpX.F90 np_force.F90  pb_fdfrc.F90 pb_force.F90 pb_init.F90
	pb_p3m.F90 pb_read.F90 pb_write.F90 sa_driver.F90 timer.F90
	curv.F90 density.F90 dsvdc.F90 gen_dx_file.F90 GrToPr.F90
	memory_module.F90 miniop.F90 pb_bldsys.F90 pb_crgview.F90
	pb_direct.F90 pb_fddrv.F90  pb_list.F90 pb_lsolver.F90  pb_mpfrc.F90
	pb_nlsolver.F90 pbsa_lib.F90 phidev.F90 project.F90 pythag.F90
	svbksb.F90 svdcmp.F90 transf.F90)

# sources for hooking up apps using libtorch
set(TORCH_SOURCES mlses.cpp utils.cpp mlses_wrapper.cpp)

# sources for libpbsa_cuda
set(PBSA_CUDA_CUDA_SOURCES mlses_predict.cu cusp_LinearSolvers.cu cusparse_LinearSolvers.cu cuda_pb.cu kLinearSolvers.cu cuda_LinearSolvers.cu)
set(PBSA_CUDA_CXX_SOURCES dummy.cpp)

#------------------------------------------------------------------------------------------

if(USE_FFT)
	add_definitions(-DFFTW)
endif()

apply_optimization_declarations(${EXECUTABLE_SOURCES} ${LIBPBSA_SOURCES} ${SANDER_PBSA_SOURCES})
set_property(SOURCE interface.c PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS})

#------------------------------------------------------------------------------------------
# copy traced mlses model to build dir
if(libtorch_ENABLED)
	file(COPY mlses_traced.pt DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
endif()

#------------------------------------------------------------------------------------------
#handle pbsa executable
if(libtorch_ENABLED)
	add_executable(pbsa ${EXECUTABLE_SOURCES} ${TORCH_SOURCES})
	target_link_libraries(pbsa netlib amber_common netcdff netcdf ${LIBTORCH_LIBRARIES})
	target_include_directories(pbsa PRIVATE ${LIBTORCH_INCLUDE_DIRS})
	target_compile_definitions(pbsa PRIVATE TORCH)
else()
	add_executable(pbsa ${EXECUTABLE_SOURCES})
	target_link_libraries(pbsa netlib amber_common netcdff netcdf)
endif()

config_module_dirs(pbsa ${PBSAEXE_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

if(USE_FFT)
	target_link_libraries(pbsa fftw)
endif()

set_target_properties(pbsa PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

#------------------------------------------------------------------------------------------
#libpbsa
add_library(libpbsa ${LIBPBSA_SOURCES})
set_target_properties(libpbsa PROPERTIES COMPILE_DEFINITIONS LIBPBSA)
config_module_dirs(libpbsa ${LIBPBSA_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
remove_prefix(libpbsa)

#add libraries to interface link libraries
if(USE_FFT)
	target_link_libraries(libpbsa fftw)
endif()

target_link_libraries(libpbsa netlib amber_common)

install_libraries(libpbsa)

resolve_cmake_library_list(UNUSED libpbsa)

#------------------------------------------------------------------------------------------
#serial/simple pbsa
add_executable(simplepbsa simplepbsa.c gopt.c)
set_property(TARGET simplepbsa PROPERTY COMPILE_OPTIONS ${OPT_CFLAGS}) #not sure whether this is supposed to be optimized... airing towards yes
target_link_libraries(simplepbsa libpbsa netlib sff amber_common)

install(TARGETS pbsa simplepbsa DESTINATION ${BINDIR})

if(libtorch_ENABLED)
	install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/mlses_traced.pt DESTINATION ${LIBDIR})
endif()

#------------------------------------------------------------------------------------------
#SANDER pbsa library

#cannot be built as shared because it has referenes to sander functions
add_library(pbsa_sander STATIC ${SANDER_PBSA_SOURCES})

set_target_properties(pbsa_sander PROPERTIES COMPILE_DEFINITIONS SANDER)

config_module_dirs(pbsa_sander ${LIBPBSA_SANDER_MOD_DIR} ${SANDER_COMMON_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
add_dependencies(pbsa_sander sander_base_obj) #force compilation of needed modules
set_property(TARGET pbsa_sander PROPERTY POSITION_INDEPENDENT_CODE ${SHARED})
target_link_libraries(pbsa_sander netlib amber_common)

install_libraries(pbsa_sander)

# --------------------------------------------------------------------
# Parallel versions
# --------------------------------------------------------------------

if(MPI)
	make_mpi_version(libpbsa libpbsa_mpi LANGUAGES Fortran INSTALL)
	config_module_dirs(libpbsa_mpi ${LIBPBSA_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

	make_mpi_version(simplepbsa simplepbsa.MPI LANGUAGES Fortran INSTALL)

	make_mpi_version(pbsa_sander pbsa_sander_mpi LANGUAGES Fortran)
	add_dependencies(pbsa_sander sander_base_obj_mpi libpbsa) #force compilation of needed modules
	config_module_dirs(pbsa_sander_mpi ${LIBPBSA_SANDER_MPI_MOD_DIR} ${LIBPBSA_MPI_MOD_DIR} ${SANDER_COMMON_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
endif()

if(CUDA)
	if(NOT HIP)
		include_directories(../cusplibrary-cuda9)
	else()
		include_directories(../cusplibrary-hip)
		include_directories(${CUDA_TOOLKIT_ROOT_DIR}/rocthrust/include)
	endif()
	
	# set cuda_cg-wrapper.cu options
	set(PBSA_CUDA_DEFINITIONS
		-DCUSP # select underlying library -> cusp [CUSPARSE, CUSP]
		-DPCG # select method -> PCG [CG, PCG]
		-DJacobi # [ILU0, IC02, Jacobi, Smooth]
		-DDIA # select matrix format -> DIA [CSR, DIA, ELL, COO, HYB]
		-DCUDA_MG) # [DCUDA_MG, NULL - for other CUDA solvers]

	if(libtorch_ENABLED AND CUDNN)
		cuda_add_executable(pbsa.cuda ${EXECUTABLE_SOURCES} ${TORCH_SOURCES} ${PBSA_CUDA_CUDA_SOURCES} OPTIONS ${PBSA_CUDA_DEFINITIONS})
		target_link_libraries(pbsa.cuda cublas cusparse netlib amber_common netcdff netcdf ${LIBTORCH_LIBRARIES})
		target_include_directories(pbsa.cuda PRIVATE ${LIBTORCH_INCLUDE_DIRS})

		# tell the Fortran code to use the cuda versions of routines
		# Matrix format macro DIA/CSR is also required to make pb_lsolver.F90 correctly compiled, DIA by default
		target_compile_definitions(pbsa.cuda PRIVATE nvCUDA DIA TORCH)
	else()
		cuda_add_executable(pbsa.cuda ${EXECUTABLE_SOURCES} ${PBSA_CUDA_CUDA_SOURCES} ${PBSA_CUDA_CXX_SOURCES} OPTIONS ${PBSA_CUDA_DEFINITIONS})
		target_link_libraries(pbsa.cuda cublas cusparse netlib amber_common netcdff netcdf)

		# tell the Fortran code to use the cuda versions of routines
		# Matrix format macro DIA/CSR is also required to make pb_lsolver.F90 correctly compiled, DIA by default
		target_compile_definitions(pbsa.cuda PRIVATE nvCUDA DIA)
	endif()

	config_module_dirs(pbsa.cuda ${PBSAEXE_CUDA_MOD_DIR} ${AMBER_COMMON_MOD_DIR})

	# imposing higher standard for wrapping libtorch
	set_target_properties(pbsa.cuda PROPERTIES CXX_STANDARD 14 CXX_STANDARD_REQUIRED YES)

	# Make sure the Fortran main gets linked
	set_property(TARGET pbsa.cuda PROPERTY LINKER_LANGUAGE Fortran)

	install(TARGETS pbsa.cuda DESTINATION ${BINDIR} COMPONENT CUDA)
endif()
