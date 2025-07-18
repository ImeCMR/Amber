#CMake buildfile for PMEMD
add_subdirectory(xray)

set(PMEMD_FORTRAN_SOURCES gbl_constants.F90 gbl_datatypes.F90 state_info.F90 file_io_dat.F90
        external_dat.F90
        mdin_ctrl_dat.F90 mdin_emil_dat.F90 mdin_ewald_dat.F90 mdin_debugf_dat.F90 prmtop_dat.F90
        inpcrd_dat.F90 dynamics_dat.F90 emil.F90 img.F90 nbips.F90 offload_allocation.F90
        parallel_dat.F90 parallel.F90 gb_parallel.F90
        pme_direct.F90 pme_recip_dat.F90 pme_slab_recip.F90 pme_blk_recip.F90
        pme_slab_fft.F90 pme_blk_fft.F90 pme_fft_dat.F90 fft1d.F90
        bspline.F90 nebread.F90 neb.F90 pme_force.F90 pbc.F90 nb_pairlist.F90 gb_ene_hybrid.F90
        nb_exclusions.F90 cit.F90 dynamics.F90 bonds.F90 angles.F90 dihedrals.F90
        extra_pnts_nb14.F90 runmd.F90 resamplekin.F90 loadbal.F90 shake.F90 prfs.F90 mol_list.F90
        runmin.F90 constraints.F90 axis_optimize.F90 gb_ene.F90 veclib.F90 gb_force.F90
        timers.F90 pmemd_lib.F90 runfiles.F90 file_io.F90
        bintraj.F90 binrestart.F90 phmd.F90
        pmemd.F90 random.F90 degcnt.F90 erfcfun.F90 nmr_calls.F90 nmr_lib.F90
        get_cmdline.F90 master_setup.F90 pme_alltasks_setup.F90 pme_setup.F90
        ene_frc_splines.F90 gb_alltasks_setup.F90 nextprmtop_section.F90
        angles_ub.F90 dihedrals_imp.F90 cmap.F90 charmm.F90 charmm_gold.F90
        findmask.F90 remd.F90 multipmemd.F90 remd_exchg.F90 amd.F90 gamd.F90 ti.F90 gbsa.F90
        barostats.F90 scaledMD.F90 constantph.F90 constante.F90 energy_records.F90 constantph_dat.F90
        constante_dat.F90 external.F90 relaxmd.F90 sgld.F90 emap.F90 get_efield_energy.F90
		processor.F90 parallel_processor.F90 pme_fft_midpoint.F90 pme_recip_midpoint.F90
        dihedrals_midpoint.F90 angles_midpoint.F90 bonds_midpoint.F90 runreweight.F90
        dbg_arrays.F90 mcres.F90 ensure_alloc.F90 reservoir.F90 sams.F90 ramd.F90
        asm_CV.F90 asm_multiCV.F90 asm_splines_utilities.F90
        asm_utilities.F90 asm_CV_utilities.F90 asm.F90
        gti.F90 reaf.F90 rmsd.F90
        xray_globals.F90
        xray_interface.F90
        xray_interface_impl_cpu.F90
        xray_interface_pre_init_data.F90
        constants.F90 assert.F90
        hybridsolvent_remd.F90
        md_scheme.F90)

set(PMEMD_FORTRAN_GPU_SOURCES
    xray_interface_impl_gpu.F90
)

set(PMEMD_C_SOURCES pmemd_clib.c)

set(PMEMD_CXX_SOURCES boost-gamma-distribution.cpp)

# PLUMED
set(PLUMED_SOURCE Plumed.c)

#Fortran NFE sources
set(PMEMD_NFE_SOURCES nfe_lib.F90 nfe_setup.F90 nfe_colvar.F90 nfe_smd.F90 nfe_abmd.F90 nfe_pmd.F90 nfe_bbmd.F90 nfe_stsm.F90)

message(STATUS "KMMD_LIB: ${KMMD_LIB}")

#build settings
#------------------------------------------------------------------------------------------

set_property(SOURCE ${PMEMD_C_SOURCES} PROPERTY COMPILE_FLAGS "${PMEMD_CFLAGS_SPC}")

set_property(SOURCE ${PMEMD_FORTRAN_SOURCES} ${PMEMD_NFE_SOURCES} PROPERTY COMPILE_FLAGS "${PMEMD_FFLAGS_SPC}")

# compile pmemd prmtop_dat at lower optimization for buggy gnu 5.x: see bug 303.
if(PMEMD_GNU_BUG_303)
	append_compile_flags(-fno-tree-vectorize prmtop_dat.F90)
endif()

# PLUMED
set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_FLAGS "${PMEMD_NO_OPT_CFLAGS_SPC}")

if(plumed_ENABLED)
    set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_DEFINITIONS __PLUMED_STATIC_KERNEL)
elseif(PLUMED_RUNTIME_LINK)
    set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_DEFINITIONS __PLUMED_HAS_DLOPEN)
endif()

#executables
#------------------------------------------------------------------------------------------

include_directories(.)

if(CUDA)
	add_subdirectory(cuda)

	# avoid building all of the source code multiple times by making an object library
	# --------------------------------------------------------------------
	add_library(pmemd_obj_cuda OBJECT ${PMEMD_C_SOURCES} ${PMEMD_CXX_SOURCES} ${PMEMD_FORTRAN_SOURCES} ${PMEMD_FORTRAN_GPU_SOURCES} ${PMEMD_NFE_SOURCES} ${PLUMED_SOURCE})
	config_module_dirs(pmemd_obj_cuda ${PMEMD_MOD_DIR}/cuda ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})

	add_dependencies(pmemd_obj_cuda amber_common netcdff)

	target_compile_definitions(pmemd_obj_cuda PRIVATE ${PMEMD_DEFINITIONS} ${PMEMD_CUDA_DEFINES})
    target_include_directories(pmemd_obj_cuda PRIVATE $<TARGET_PROPERTY:boost_headers,INTERFACE_INCLUDE_DIRECTORIES>)

    if(${CMAKE_VERSION} VERSION_LESS "3.12.0")
        target_include_directories(
                pmemd_obj_cuda
                   PUBLIC
                $<TARGET_PROPERTY:PMEMD::xray_gpu,INTERFACE_INCLUDE_DIRECTORIES>
        )
        target_compile_definitions(
                pmemd_obj_cuda
                    PUBLIC
                $<TARGET_PROPERTY:PMEMD::xray_gpu,INTERFACE_COMPILE_DEFINITIONS>
        )
        add_dependencies(pmemd_obj_cuda PMEMD::xray_gpu)
    else()
        target_link_libraries(pmemd_obj_cuda PUBLIC PMEMD::xray_gpu) # requires CMake-3.12
    endif()

	if(MPI)
		make_mpi_version(pmemd_obj_cuda pmemd_obj_cuda_mpi LANGUAGES Fortran)
        target_compile_options(pmemd_obj_cuda_mpi PRIVATE ${PMEMD_MPI_FLAGS})
        target_compile_definitions(pmemd_obj_cuda_mpi PRIVATE ${PMEMD_MPI_DEFINITIONS})
		config_module_dirs(pmemd_obj_cuda_mpi ${PMEMD_MOD_DIR}/cuda_mpi ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})
	endif()

	# now make the executables
	# --------------------------------------------------------------------
	foreach(PRECISION ${PMEMD_CUDA_PRECISIONS})

		if(HIP)
			set(EXE_NAME pmemd.hip_${PRECISION})
		else()
			set(EXE_NAME pmemd.cuda_${PRECISION})
		endif()

		add_executable(${EXE_NAME} $<TARGET_OBJECTS:pmemd_obj_cuda>)
		set_property(TARGET ${EXE_NAME} PROPERTY LINKER_LANGUAGE Fortran)
        target_link_libraries(${EXE_NAME} pmemd_cuda_${PRECISION} emil amber_common netcdff netlib boost_headers ${KMMD_LIB} PMEMD::xray_gpu)
		cuda_add_cufft_to_target(${EXE_NAME})


        # PLUMED
        if(PLUMED_RUNTIME_LINK)
            target_link_libraries(${EXE_NAME} dl)
            set_target_properties(${EXE_NAME} PROPERTIES ENABLE_EXPORTS TRUE)
        else()
            if(plumed_ENABLED)
                target_link_libraries(${EXE_NAME} plumed)
            endif()
        endif()

		install(TARGETS ${EXE_NAME} DESTINATION ${BINDIR} COMPONENT pmemd_CUDA)

		if(MPI)
			add_executable(${EXE_NAME}.MPI $<TARGET_OBJECTS:pmemd_obj_cuda_mpi>)
			set_property(TARGET ${EXE_NAME}.MPI PROPERTY LINKER_LANGUAGE Fortran)
			target_link_libraries(${EXE_NAME}.MPI pmemd_cuda_${PRECISION}_mpi emil_mpi amber_common netcdff netlib mpi_c mpi_cxx mpi_fortran boost_headers PMEMD::xray_gpu ${KMMD_LIB})

            # PLUMED
            if(PLUMED_RUNTIME_LINK)
                target_link_libraries(${EXE_NAME}.MPI dl)
                set_target_properties(${EXE_NAME}.MPI PROPERTIES ENABLE_EXPORTS TRUE)
            else()
                if(plumed_ENABLED)
                    target_link_libraries(${EXE_NAME}.MPI plumed)
                endif()
            endif()

			install(TARGETS ${EXE_NAME}.MPI DESTINATION ${BINDIR} COMPONENT pmemd_CUDA)
		endif()

	endforeach()

	if(HIP)
	#copy pmemd.hip.SPFP to pmemd.hip.
	install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.hip_${PMEMD_DEFAULT_PRECISION}
		\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.hip)" COMPONENT pmemd_CUDA)

        if(MPI)
                install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        \$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.hip_${PMEMD_DEFAULT_PRECISION}.MPI
			\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.hip.MPI)" COMPONENT pmemd_CUDA)
        endif()

        else()
	#copy pmemd.cuda.SPFP to pmemd.cuda.
	#This is actually kind of tricky, as CMake doesn't support renaming targets on install, so we have to do it ourselves.
	install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
		\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.cuda_${PMEMD_DEFAULT_PRECISION}
		\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.cuda)" COMPONENT pmemd_CUDA)

	if(MPI)
		install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
			\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.cuda_${PMEMD_DEFAULT_PRECISION}.MPI
			\$ENV{DESTDIR}\${CMAKE_INSTALL_PREFIX}/${BINDIR}/pmemd.cuda.MPI)" COMPONENT pmemd_CUDA)
	endif()
	endif()
endif()

add_executable(pmemd ${PMEMD_C_SOURCES} ${PMEMD_CXX_SOURCES} ${PMEMD_FORTRAN_SOURCES} ${PMEMD_NFE_SOURCES} ${PLUMED_SOURCE})
config_module_dirs(pmemd ${PMEMD_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})

# Address sanitizer must be first linked library
target_link_libraries(pmemd emil amber_common netcdff netlib boost_headers ${KMMD_LIB} PMEMD::xray_cpu)
target_compile_definitions(pmemd PRIVATE ${PMEMD_DEFINITIONS})
set_property(TARGET pmemd PROPERTY LINKER_LANGUAGE Fortran)


# PLUMED
if(PLUMED_RUNTIME_LINK)
    target_link_libraries(pmemd dl)
    set_target_properties(pmemd PROPERTIES ENABLE_EXPORTS TRUE)
else()
    if(plumed_ENABLED)
        target_link_libraries(pmemd plumed::plumed)
    endif()
endif()

if(mbx_ENABLED)
    target_link_libraries(pmemd MBX::mbx)
    target_compile_definitions(pmemd PRIVATE MBX)
endif()

install(TARGETS pmemd DESTINATION ${BINDIR} COMPONENT pmemd)

# MPI parallelization
if(MPI)
	make_mpi_version(pmemd pmemd.MPI LANGUAGES Fortran)
	config_module_dirs(pmemd.MPI ${PMEMD_MOD_DIR}/mpi ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})

    target_compile_definitions(pmemd.MPI PRIVATE ${PMEMD_MPI_DEFINITIONS})
    target_compile_options(pmemd.MPI PRIVATE ${PMEMD_MPI_FLAGS})

	install(TARGETS pmemd.MPI DESTINATION ${BINDIR} COMPONENT pmemd_MPI)

    if(OPENMP AND PMEMD_OMP_MPI)
        make_openmp_version(pmemd.MPI pmemd.OMP.MPI LANGUAGES Fortran)
        config_module_dirs(pmemd.OMP.MPI ${PMEMD_MOD_DIR}/omp_mpi ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})
        
        target_compile_definitions(pmemd.OMP.MPI PRIVATE ${PMEMD_OMP_DEFINITIONS})
        target_compile_options(pmemd.OMP.MPI PRIVATE ${PMEMD_OMP_FLAGS})
        target_link_libraries(pmemd.OMP.MPI openmp_fortran ${KMMD_LIB})

        install(TARGETS pmemd.OMP.MPI DESTINATION ${BINDIR} COMPONENT pmemd_MPI)
    endif()

endif()

# MIC parallelization
if(MIC_KL)

    # mic2 only exists as an MPI version
    copy_target(pmemd.MPI pmemd.mic_kl.MPI)
    config_module_dirs(pmemd.mic_kl.MPI ${PMEMD_MOD_DIR}/mic_kl ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR})
    target_link_libraries(pmemd.mic_kl.MPI openmp_fortran)

    # add mic flags
    target_compile_options(pmemd.mic_kl.MPI PRIVATE ${PMEMD_MIC_FLAGS})

    # also add openmp flags
    target_compile_definitions(pmemd.mic_kl.MPI PRIVATE ${PMEMD_OMP_DEFINITIONS})
    target_compile_options(pmemd.mic_kl.MPI PRIVATE ${PMEMD_OMP_FLAGS})

    install(TARGETS pmemd.mic_kl.MPI DESTINATION ${BINDIR} COMPONENT pmemd_MIC)
endif()



