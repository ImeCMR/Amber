#AMBER CMake buildfile for sander

#Source Listing
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#common source.
#Contained in every executable and library, but may get built differently in different versions
set(SANDER_COMMON_SOURCE resamplekin.F90 boost-gamma-distribution.cpp les.F90 state.F90  cmd_vars.F90 nose_hoover.F90
	abfqmmm_module.F90 charmm.F90 molecule.F90 parms.F90 memory_module.F90 stack.F90
    trace.F90 rmsgrd.F90  icosasurf.F90 findmask.F90  ips.F90 file_io_dat.F90 decomp.F90
    softcore.F90 ew_bspline.F90 bintraj.F90 sglds.F90  nonbond_list.F90
    amoeba_valence.F90 amoeba_multipoles.F90 amoeba_recip.F90 amoeba_interface.F90
    amoeba_direct.F90 amoeba_mdin.F90 amoeba_adjust.F90 amoeba_self.F90
    amoeba_vdw.F90 amoeba_induced.F90 amoeba_runmd.F90
    pGM_multipoles.F90 pGM_recip.F90 pGM_interface.F90
    pGM_direct.F90 pGM_mdin.F90 pGM_self.F90 pGM_induced.F90
    pimd_vars.F90 pimd_force.F90 pimd_init.F90 cmd_matrix.F90
    lscivr_vars.F90 lsc_init.F90 lsc_xp.F90
    nfe-abmd-hooks.F90 nfe-bbmd-ctxt.F90 nfe-bbmd-hooks.F90
    nfe-colvar.F90 nfe-colvar-math.F90 nfe-colvar-type.F90 nfe-colvar-utils.F90
    nfe-constants.F90 nfe-cv-ANGLE.F90 nfe-cv-COM_ANGLE.F90
    nfe-cv-COM_DISTANCE.F90 nfe-cv-COM_TORSION.F90 nfe-cv-COS_OF_DIHEDRAL.F90
    nfe-cv-DISTANCE.F90 nfe-cv-HANDEDNESS.F90 nfe-cv-LCOD.F90
    nfe-cv-MULTI_RMSD.F90 nfe-cv-N_OF_BONDS.F90 nfe-cv-N_OF_STRUCTURES.F90
    nfe-cv-R_OF_GYRATION.F90 nfe-cv-TORSION.F90
    nfe-cv-ORIENTATION_ANGLE.F90 nfe-cv-ORIENTATION_PROJ.F90
    nfe-cv-SPINANGLE.F90 nfe-cv-TILT.F90
    nfe-cv-QUATERNION0.F90 nfe-cv-QUATERNION1.F90
    nfe-cv-QUATERNION2.F90 nfe-cv-QUATERNION3.F90
    nfe-pmd-hooks.F90 nfe-rmsd.F90 nfe-sander-hooks.F90 nfe-sander-proxy.F90
    nfe-stsm-hooks.F90 nfe-cv-PCA.F90 nfe-read-pca.F90 nfe-cv-SIN_OF_DIHEDRAL.F90
    nfe-cv-PAIR_DIHEDRAL.F90 nfe-cv-PATTERN_DIHEDRAL.F90 nfe-cv-DF_COM_DISTANCE.F90
    nfe-smd-hooks.F90 nfe-umbrella.F90 nfe-utils.F90 nfe-value.F90 mt19937.F90
    pol_gauss.F90 quench.F90 external.F90 relaxmd.F90 music_module.F90 constante.F90
    dprcmod.F90 dprc.cpp f90stream.cpp f90streammod.F90
    asm_CV.F90 asm_multiCV.F90 asm_splines_utilities.F90 asm_utilities.F90 asm_CV_utilities.F90 asm.F90
    rattle.F90 md_scheme.F90 cewmod.F90)

set(EMIL_SOURCE emil.F90 mdin_emil_dat.F90)

set(SEBOMD_SOURCE sebomd_module.F90 sebomd_arrays.F90)

#Sources only built into certain executables
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

set(EVB_BASE_SOURCE evb_vars.F90 evb_mathf.F90 evb_init.F90 evb_bcast.F90
    evb_dealloc.F90 check_input.F90 evb_io.F90 evb_keywrd.F90
    exchange_gauss.F90 exchange_warshel.F90 react_flux_init.F90
    react_path.F90 react_flux.F90 kappa_keck.F90
    schlegel_dg.F90 schlegel_full.F90 schlegel_diis.F90 schlegel_gmres.F90
    schlegel_gbasis.F90 schlegel_pgbasis.F90
    bnd_strtch.F90 ang_bend.F90 dihed_tors.F90 torinternal.F90
    evb_gfchk.F90 bmat_inv.F90 bmat_grad.F90
    unitv.F90 crossv.F90 cart2internal.F90 wdc_bmat.F90 select_coord.F90
    schlegel_uff.F90 schlegel_dihed.F90 schlegel_angle.F90
    schlegel_poh_uff.F90 schlegel_poh_irc.F90 schlegel_hcn_cart.F90
    schlegel_hcn_int_irc.F90 schlegel_hcn_cart_irc.F90 evb_outoptions.F90)

set(EVB_DEBUG_SOURCE evb_2stdebug.F90 egap_umb_2stdebug.F90 dg_debug.F90
    evb_io_debug.F90 bmat_debug.F90 bmat_inv_debug.F90 evb_io_debug2.F90
    morse_anal2num.F90 xwarshel_anal2num.F90 xgauss_anal2num.F90
    vdw_anal2num.F90 rcdbonds_anal2num.F90 rcbond_anal2num.F90
    dg_grad_anal2num.F90)

set(EVB_SOURCE ${EVB_BASE_SOURCE} ${EVB_DEBUG_SOURCE} evb_input.F90 evb_alloc.F90 evb_matrix.F90
    morsify.F90 mod_vdw.F90 evb_ntrfc.F90 evb_force.F90 out_evb.F90
    dbonds_anal2num.F90 bond_anal2num.F90 evb_umb.F90 evb_umb_primitive.F90)

set(EVB_PIMD_SOURCE ${EVB_BASE_SOURCE} ${EVB_DEBUG_SOURCE} evb_input.F90 evb_alloc.F90 evb_matrix.F90
    morsify.F90 mod_vdw.F90 evb_ntrfc.F90 evb_force.F90
    evb_pimd_init.F90 out_evb.F90 dbonds_anal2num.F90 bond_anal2num.F90
    evb_umb.F90 evb_umb_primitive.F90)

#main source for sander.APBS
set(APBS_SOURCE nose_hoover_init.F90
    lmod.F90 egb.F90 apbs_vars.F90 apbs.F90
    remd.F90 force.F90
    relax_mat.F90 nmr.F90 multitmd.F90
    multisander.F90 sander.F90 trajene.F90 cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90
    printe.F90 runmin.F90 rdparm.F90
    mdread.F90 locmem.F90 runmd.F90 getcor.F90 degcnt.F90 decnvh.F90
    fastwt.F90 parallel.F90 shake.F90 ene.F90 mdwrit.F90 minrit.F90
    set.F90 dynlib.F90 mdfil.F90 nmlsrc.F90
    ew_force.F90 ew_setup.F90 ew_box.F90
    ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90 align.F90 csa.F90
    rfree.F90 rgroup.F90 amopen.F90
    debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
    mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90
    matinv.F90 assert.F90 erfcfun.F90 sinr_t.F90
    constantph.F90 prn_dipoles.F90  iwrap2.F90 emap.F90
    binrestart.F90
    spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90 dssp.F90
    crg_reloc.F90 sander_lib.F90 amd.F90 scaledMD.F90 barostats.F90
    fsockets_pure.F90 ipimod.F90)

set(QM_SOURCE qm_mm.F90 qm_ewald.F90 qmmm_vsolv.F90
    qm2_extern_util_module.F90 qm2_extern_module.F90 qm2_extern_adf_module.F90 qm2_read_adf_results.c KFReader.c
    ArrayList.c qm2_extern_gms_module.F90 qm2_extern_tc_module.F90 qm2_extern_gau_module.F90 qm2_extern_mrcc_module.F90
    qm2_extern_orc_module.F90 qm2_extern_nw_module.F90 qm2_extern_qc_module.F90 qm2_extern_reaxff_puremd_module.F90
    qm2_extern_quick_module.F90 qm2_extern_genmpi_module.F90
    qm_ewald_comm.F90 xtb_module.F90 dftbplus_module.F90
    qmmm_adaptive_module.F90 prn_qmmm_dipole.F90 quick_module.F90 tcpb_module.F90)

set(QM_API_SOURCE qm_ewald.F90 qmmm_vsolv.F90
    qm2_extern_util_module.F90 qm2_extern_module.F90 qm2_extern_adf_module.F90 qm2_read_adf_results.c KFReader.c
    ArrayList.c qm2_extern_gms_module.F90 qm2_extern_tc_module.F90 qm2_extern_gau_module.F90 qm2_extern_mrcc_module.F90
    qm2_extern_orc_module.F90 qm2_extern_nw_module.F90 qm2_extern_qc_module.F90 qm2_extern_reaxff_puremd_module.F90
    qm2_extern_quick_module.F90 qm2_extern_genmpi_module.F90
    qm_ewald_comm.F90 xtb_module.F90 dftbplus_module.F90
    qmmm_adaptive_module.F90 prn_qmmm_dipole.F90 quick_module.F90 tcpb_module.F90)

set(API_SOURCE  nose_hoover_init.F90
    lmod.F90 egb.F90
    remd.F90 force.F90
    relax_mat.F90 nmr.F90 multitmd.F90
    trajene.F90 cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90
    printe.F90 runmin.F90 rdparm.F90 cpparm.F90 sander_lib.F90
    mdread.F90 locmem.F90 runmd.F90 getcor.F90 degcnt.F90 decnvh.F90
    fastwt.F90 parallel.F90 shake.F90 ene.F90 mdwrit.F90 minrit.F90
    set.F90 dynlib.F90 mdfil.F90 nmlsrc.F90
    ew_force.F90 ew_setup.F90 ew_box.F90
    ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90 align.F90 csa.F90
    rfree.F90 rgroup.F90 amopen.F90
    debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
    mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90
    matinv.F90 assert.F90 erfcfun.F90 sinr_t.F90
    constantph.F90 prn_dipoles.F90  iwrap2.F90 emap.F90 binrestart.F90
    spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90 dssp.F90
    crg_reloc.F90 amd.F90 scaledMD.F90 barostats.F90 interface.F90  qm_mm.F90 prmtop_type.F90
    fsockets_pure.F90 ipimod.F90)

set(LES_API_SOURCE egb.F90 nose_hoover_init_LES.F90
    remd.F90 relax_mat.F90 nmr.F90 multitmd.F90
    cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90 cpparm.F90
    printe.F90 runmin.F90 force.F90 rdparm.F90
    mdread.F90 locmem.F90 runmd.F90 degcnt.F90 decnvh.F90
    fastwt.F90 parallel.F90 shake.F90 ene.F90
    dynlib.F90 nmlsrc.F90
    ew_force.F90 ew_setup.F90
    ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90
    align.F90 csa.F90 rfree.F90 rgroup.F90
    amopen.F90 sander_lib.F90
    debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
    mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90 matinv.F90 assert.F90
    erfcfun.F90 sinr_t.F90
    constantph.F90 prn_dipoles.F90 iwrap2.F90 emap.F90
    spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90 dssp.F90
    crg_reloc.F90 barostats.F90 prmtop_type.F90
    binrestart.F90 getcor.F90 mdfil.F90 mdwrit.F90
    lmod.F90 trajene.F90 ew_box.F90 minrit.F90 set.F90 qm_mm.F90
    amd.F90 scaledMD.F90 interface.F90
    fsockets_pure.F90 ipimod.F90)


#main source for sander.PUPIL
set(PUPIL_SOURCE nose_hoover_init.F90
	lmod.F90 egb.F90
	pupildata.F90
	relax_mat.F90 nmr.F90 multitmd.F90
	multisander.F90 sander.F90 trajene.F90 cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90
	printe.F90 runmin.F90 force.F90 rdparm.F90
	mdread.F90 locmem.F90 runmd.F90 getcor.F90 degcnt.F90 decnvh.F90
	fastwt.F90 parallel.F90 shake.F90 ene.F90 mdwrit.F90 minrit.F90
	set.F90 dynlib.F90 mdfil.F90 nmlsrc.F90
	ew_force.F90 ew_setup.F90 ew_box.F90
	ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90 align.F90 csa.F90
	rfree.F90 rgroup.F90 amopen.F90 sander_lib.F90
	debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
	mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90
	matinv.F90 assert.F90 erfcfun.F90 sinr_t.F90
	constantph.F90 prn_dipoles.F90  iwrap2.F90 emap.F90 binrestart.F90
	spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90
	putvalues.F90 fixport.F90  dssp.F90 crg_reloc.F90 amd.F90 scaledMD.F90
	barostats.F90
	fsockets_pure.F90 ipimod.F90)

set(LES_SOURCE nose_hoover_init_LES.F90
    remd.F90 egb.F90
    relax_mat.F90 nmr.F90 multitmd.F90
    multisander.F90 sander.F90  cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90
    printe.F90 runmin.F90 force.F90 rdparm.F90
    mdread.F90 locmem.F90 runmd.F90 degcnt.F90 decnvh.F90
    fastwt.F90 parallel.F90 shake.F90 ene.F90
    dynlib.F90 nmlsrc.F90
    ew_force.F90 ew_setup.F90
    ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90
    align.F90 csa.F90 rfree.F90 rgroup.F90
    amopen.F90 sander_lib.F90
    debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
    mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90 matinv.F90 assert.F90
    erfcfun.F90 sinr_t.F90
    constantph.F90 prn_dipoles.F90  iwrap2.F90 emap.F90
    spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90 dssp.F90
    crg_reloc.F90 barostats.F90
    qm_ewald.F90 qmmm_vsolv.F90
    qm2_extern_util_module.F90 qm2_extern_module.F90 qm2_extern_adf_module.F90 qm2_read_adf_results.c KFReader.c
    ArrayList.c qm2_extern_gms_module.F90 qm2_extern_tc_module.F90 qm2_extern_gau_module.F90
    qm2_extern_orc_module.F90 qm2_extern_nw_module.F90 qm2_extern_mrcc_module.F90
    qm2_extern_qc_module.F90 qm2_extern_reaxff_puremd_module.F90
    qm_ewald_comm.F90 xtb_module.F90 dftbplus_module.F90
    qm2_extern_quick_module.F90 qm2_extern_genmpi_module.F90 qmmm_adaptive_module.F90
    prn_qmmm_dipole.F90 binrestart.F90 getcor.F90 mdfil.F90 mdwrit.F90
    lmod.F90 trajene.F90 ew_box.F90 minrit.F90 set.F90 qm_mm.F90 amd.F90 scaledMD.F90 quick_module.F90 tcpb_module.F90
    fsockets_pure.F90 ipimod.F90)

#mm code and miscellaneous files
set(MM_SOURCE  force.F90 nose_hoover_init.F90
    lmod.F90 egb.F90
    remd.F90
    relax_mat.F90 nmr.F90 multitmd.F90
    multisander.F90 sander.F90 trajene.F90 cshf.F90 nmrcal.F90 mtmdcall.F90 pearsn.F90
    printe.F90 runmin.F90 rdparm.F90 sander_lib.F90
    mdread.F90 locmem.F90 runmd.F90 getcor.F90 degcnt.F90 decnvh.F90
    fastwt.F90 parallel.F90 shake.F90 ene.F90 mdwrit.F90 minrit.F90
    set.F90 dynlib.F90 mdfil.F90 nmlsrc.F90
    ew_force.F90 ew_setup.F90 ew_box.F90
    ew_fft.F90 short_ene.F90 ew_recip.F90 pcshift.F90 align.F90 csa.F90
    rfree.F90 rgroup.F90 amopen.F90
    debug.F90 ew_recip_reg.F90 ew_handle_dips.F90 ew_bndpol.F90 ew_dipole_recip.F90
    mexit.F90 new_time.F90 extra_pts.F90 thermo_int.F90
    matinv.F90 assert.F90 erfcfun.F90 sinr_t.F90
    constantph.F90 prn_dipoles.F90  iwrap2.F90 emap.F90 binrestart.F90
    spatial_recip.F90 spatial_fft.F90  mbar.F90 linear_response.F90 xref.F90 dssp.F90
    crg_reloc.F90 amd.F90 scaledMD.F90 barostats.F90
    fsockets_pure.F90 ipimod.F90)

set(AMBMASK_SOURCE ambmask.F90 findmask.F90 amopen.F90 parms.F90 memory_module.F90 mexit.F90)

#C source for plumed
set(PLUMED_SOURCE Plumed.c)

if(USE_FFT)
	#tell sander to enable rism interop
	add_definitions(-DRISMSANDER)

	set(RISM_SANDER_IFACE $<TARGET_OBJECTS:sander_rism_interface>)
	set(RISM_SANDER_IFACE_MPI $<TARGET_OBJECTS:sander_rism_interface_mpi>)
else()
	set(RISM_SANDER_IFACE "")
	set(RISM_SANDER_IFACE_MPI "")
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#List of all Fortran sources, for setting flags
set(ALL_FORTRAN_SOURCES ${SANDER_COMMON_SOURCE} ${LES_SOURCE} ${PUPIL_SOURCE} ${LES_API_SOURCE}
	${API_SOURCE} ${MM_SOURCE} ${QM_SOURCE} ${QM_API_SOURCE} ${APBS_SOURCE} ${EVB_SOURCE} ${EVB_PIMD_SOURCE}
	${SEBOMD_SOURCE} ${AMBMASK_SOURCE} ${EMIL_SOURCE})

apply_optimization_declarations(${ALL_FORTRAN_SOURCES})

# reaxff-puremd QMMM extern interface module requires F2003 standard or later
set_source_files_properties(qm2_extern_reaxff_puremd_module.F90 PROPERTIES LANGUAGE Fortran COMPILE_FLAGS -std=f2003)

set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_FLAGS ${OPT_CFLAGS_SPC})

if(plumed_ENABLED)
	set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_DEFINITIONS __PLUMED_STATIC_KERNEL)
elseif(PLUMED_RUNTIME_LINK)
	set_property(SOURCE ${PLUMED_SOURCE} PROPERTY COMPILE_DEFINITIONS __PLUMED_HAS_DLOPEN)
endif()

if(lio_ENABLED)
    list(APPEND LES_SOURCE qm2_extern_lio_module.F90)
    list(APPEND QM_SOURCE qm2_extern_lio_module.F90)
    list(APPEND QM_API_SOURCE qm2_extern_lio_module.F90)
endif()

# EMIL is always enabled, as far as I know
add_definitions(-DEMIL)

if(mbx_ENABLED)
    add_definitions(-DMBX)
endif()

if(lio_ENABLED)
    add_definitions(-DLIO)
endif()

if(BUILD_QUICK)
    # Enable QUICK for everything
    add_definitions(-DQUICK)
    add_definitions(-DCEW)
endif()

if(BUILD_TCPB)
    # Enable TCPB for everything
    add_definitions(-DTCPB)
endif()

if(BUILD_REAXFF_PUREMD)
    # Enable reaxff-puremd interface routines
    add_definitions(-DHAVE_REAXFF_PUREMD)
endif()


#deal with modules
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

set(APBS_MOD_DIR ${SANDER_COMMON_MOD_DIR}/apbs)
set(PUPIL_MOD_DIR ${SANDER_COMMON_MOD_DIR}/pupil)
set(LES_MOD_DIR ${SANDER_COMMON_MOD_DIR}/les)
set(MPI_MOD_DIR ${SANDER_COMMON_MOD_DIR}/mpi)
set(API_MOD_DIR ${SANDER_COMMON_MOD_DIR}/api)
set(LESAPI_MOD_DIR ${SANDER_COMMON_MOD_DIR}/lesapi)
set(AMBMASK_MOD_DIR ${SANDER_COMMON_MOD_DIR}/ambmask)
file(MAKE_DIRECTORY ${APBS_MOD_DIR} ${PUPIL_MOD_DIR} ${LES_MOD_DIR} ${API_MOD_DIR} ${LESAPI_MOD_DIR} ${AMBMASK_MOD_DIR})

include_directories(${SANDER_COMMON_MOD_DIR} )

# quick mod dirs, or empty if none
if(HIP)
    set(QUICK_GPU_TARGET_SUFFIX "hip")
else()
    set(QUICK_GPU_TARGET_SUFFIX "cuda")
endif()

if(BUILD_QUICK)
    set(QUICK_SERIAL_MOD_DIR quick/serial)
    set(QUICK_MPI_MOD_DIR quick/mpi)
    set(QUICK_CUDA_MOD_DIR quick/${QUICK_GPU_TARGET_SUFFIX})
    set(QUICK_MPI_CUDA_MOD_DIR quick/mpi_${QUICK_GPU_TARGET_SUFFIX})
else()
    set(QUICK_SERIAL_MOD_DIR )
    set(QUICK_MPI_MOD_DIR )
    set(QUICK_CUDA_MOD_DIR )
    set(QUICK_MPI_CUDA_MOD_DIR )
endif()

#add final executables and libraries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

#we add all of the generated executables to the EXECUTABLES list so that we can set properties on them

#I know that this is horribly inefficient and slow.  For each executable, it compiles the whole source again.
#Make works at the object level, so the old build system could be configured to use only certain special objects built with different definitions in an executable.
#CMake, though, works at the target level, and has no idea that the same source file, built in the same way for two different targets, generates the same object file.
#I tried to extract a set of common sources which could be compiled only once.
#However, so many source files come out differently depending on the definitions they are built with.
#Alas, the module dependencies were a tangled web of pain, and I had to change it back
#Just know that I did try to fix it, and that it can't be done without moving some of the configuration from compile time to runtime
#</rant>

set(EXECUTABLES "")
set(MPI_EXECUTABLES "")
set(CUDA_MPI_EXECUTABLES "")
set(OPENMP_EXECUTABLES "")
set(CUDA_EXECUTABLES "")
set(LIBRARIES "")

# --------------------------------------------------------------------
# serial sander

#note: some of the link libraries sander depends on use its modules.
#so, we build all of the regular sander source into an object library, go build those libraries, then link everything.

add_library(sander_base_obj OBJECT ${SANDER_COMMON_SOURCE} ${MM_SOURCE} ${QM_SOURCE} ${EMIL_SOURCE} ${SEBOMD_SOURCE} ${PLUMED_SOURCE})
config_module_dirs(sander_base_obj ${SANDER_COMMON_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
add_dependencies(sander_base_obj libsqm sff_fortran netcdff libpbsa)
target_include_directories(sander_base_obj PUBLIC $<TARGET_PROPERTY:boost_headers,INTERFACE_INCLUDE_DIRECTORIES>)

if(USE_FFT)
	add_dependencies(sander_base_obj sander_rism_interface)
endif()

add_executable(sander $<TARGET_OBJECTS:sander_base_obj> ${RISM_SANDER_IFACE})
list(APPEND EXECUTABLES sander)

# --------------------------------------------------------------------
# sander PUPIL

if(BUILD_SANDER_PUPIL)
	#sander pupil executable
	add_executable(sander.PUPIL ${SANDER_COMMON_SOURCE} ${RISM_SANDER_IFACE}
		${PUPIL_SOURCE} ${QM_SOURCE} ${EMIL_SOURCE} ${SEBOMD_SOURCE} ${PLUMED_SOURCE})
	target_link_libraries(sander.PUPIL pupil::pupil)
	set_property(TARGET sander.PUPIL APPEND PROPERTY COMPILE_DEFINITIONS PUPIL_SUPPORT)
	config_module_dirs(sander.PUPIL ${PUPIL_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	list(APPEND EXECUTABLES sander.PUPIL)
endif()

# --------------------------------------------------------------------
# sander LES

if(BUILD_SANDER_LES)
	add_executable(sander.LES ${SANDER_COMMON_SOURCE} ${RISM_SANDER_IFACE}
		${LES_SOURCE} ${EMIL_SOURCE} ${SEBOMD_SOURCE} ${PLUMED_SOURCE})
	target_compile_definitions(sander.LES PRIVATE LES)
	config_module_dirs(sander.LES ${LES_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	list(APPEND EXECUTABLES sander.LES)

endif()
# --------------------------------------------------------------------
# sander APBS

if(BUILD_SANDER_APBS)
	add_executable(sander.APBS ${SANDER_COMMON_SOURCE} ${RISM_SANDER_IFACE}
		${APBS_SOURCE} ${QM_SOURCE} ${SEBOMD_SOURCE} ${EMIL_SOURCE} ${PLUMED_SOURCE})
	target_link_libraries(sander.APBS apbs)
	target_compile_definitions(sander.APBS PRIVATE APBS)
	config_module_dirs(sander.APBS ${APBS_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	list(APPEND EXECUTABLES sander.APBS)
endif()


# --------------------------------------------------------------------
# ambmask

add_executable(ambmask ${AMBMASK_SOURCE})
config_module_dirs(ambmask ${AMBMASK_MOD_DIR} ${AMBER_COMMON_MOD_DIR})
target_link_libraries(ambmask amber_common)
install(TARGETS ambmask DESTINATION ${BINDIR})

# --------------------------------------------------------------------
# sander API
if(BUILD_SANDER_API)
	add_library(libsander ${SANDER_COMMON_SOURCE} ${RISM_SANDER_IFACE}
		${API_SOURCE} ${QM_API_SOURCE} ${EMIL_SOURCE} ${SEBOMD_SOURCE} ${PLUMED_SOURCE})
	list(APPEND LIBRARIES libsander)
	target_compile_definitions(libsander PRIVATE API)
	config_module_dirs(libsander ${API_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	target_include_directories(libsander INTERFACE .)
	remove_prefix(libsander)

	#install api modules
	install(FILES ${API_MOD_DIR}/sander_api.mod sander.h DESTINATION ${INCDIR})
endif()

# --------------------------------------------------------------------
# sander LES API

if(BUILD_SANDER_LES AND BUILD_SANDER_API)

	add_library(libsanderles ${SANDER_COMMON_SOURCE} ${RISM_SANDER_IFACE}
		${LES_API_SOURCE} ${QM_API_SOURCE} ${EMIL_SOURCE} ${PLUMED_SOURCE} ${SEBOMD_SOURCE})
	list(APPEND LIBRARIES libsanderles)
	target_compile_definitions(libsanderles PRIVATE LES API)
	config_module_dirs(libsanderles ${LESAPI_MOD_DIR} ${SQMLIB_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${LIBPBSA_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	target_include_directories(libsanderles INTERFACE .)
	remove_prefix(libsanderles)

	#install api modules
	install(FILES ${LESAPI_MOD_DIR}/sanderles_api.mod DESTINATION ${INCDIR})

endif()

# --------------------------------------------------------------------
# Parallel versions
# --------------------------------------------------------------------

if(MPI)
	make_mpi_version(sander_base_obj sander_base_obj_mpi LANGUAGES Fortran SWAP_SOURCES TO ${EVB_SOURCE})
	config_module_dirs(sander_base_obj_mpi ${SANDER_COMMON_MPI_MOD_DIR} ${SQMLIB_MPI_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_MPI_DIR} ${RISMLIB_MPI_MOD_DIR} ${LIBPBSA_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_MPI_MOD_DIR})
	add_dependencies(sander_base_obj_mpi libsqm_mpi libpbsa_mpi)

	if(USE_FFT)
		add_dependencies(sander_base_obj_mpi sander_rism_interface_mpi)
	endif()

	#sander MPI executable
	add_executable(sander.MPI $<TARGET_OBJECTS:sander_base_obj_mpi> ${RISM_SANDER_IFACE_MPI})
	target_link_libraries(sander.MPI mpi_fortran)

	list(APPEND MPI_EXECUTABLES sander.MPI)

	#sander LES MPI executable
	if(BUILD_SANDER_LES)
		make_mpi_version(sander.LES sander.LES.MPI LANGUAGES Fortran SWAP_SOURCES ${RISM_SANDER_IFACE} TO ${RISM_SANDER_IFACE_MPI} ${EVB_PIMD_SOURCE})
		config_module_dirs(sander.LES.MPI ${SANDER_COMMON_MPI_MOD_DIR}/les ${SQMLIB_MPI_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_MPI_DIR} ${RISMLIB_MPI_MOD_DIR} ${LIBPBSA_MPI_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_MPI_MOD_DIR})
		list(APPEND MPI_EXECUTABLES sander.LES.MPI)
	endif()
endif()

if(OPENMP)
	make_openmp_version(sander_base_obj sander_base_obj_omp LANGUAGES Fortran)
	config_module_dirs(sander_base_obj_omp ${SANDER_COMMON_OMP_MOD_DIR} ${SQMLIB_OMP_MOD_DIR} ${RISMLIB_SANDER_INTERFACE_DIR} ${RISMLIB_MOD_DIR} ${LIBPBSA_MOD_DIR} ${AMBER_COMMON_MOD_DIR} ${NETCDF_FORTRAN_MOD_DIR} ${QUICK_SERIAL_MOD_DIR})
	add_dependencies(sander_base_obj_omp libsqm_omp sff_fortran netcdff libpbsa)

	if(USE_FFT)
		add_dependencies(sander_base_obj_omp sander_rism_interface)
	endif()

	#sander MPI executable
	add_executable(sander.OMP $<TARGET_OBJECTS:sander_base_obj_omp> ${RISM_SANDER_IFACE})
	target_link_libraries(sander.OMP openmp_fortran)
	list(APPEND OPENMP_EXECUTABLES sander.OMP)

	targets_link_libraries(${OPENMP_EXECUTABLES} LIBRARIES libsqm_omp)

endif()

if(CUDA OR HIP)

    if(BUILD_QUICK)
        # QUICK CUDA sanders are same as regular sanders but link to CUDA QUICK
        add_executable(sander.quick.${QUICK_GPU_TARGET_SUFFIX} $<TARGET_OBJECTS:sander_base_obj> ${RISM_SANDER_IFACE})
        list(APPEND CUDA_EXECUTABLES sander.quick.${QUICK_GPU_TARGET_SUFFIX})

        if(MPI)
            add_executable(sander.quick.${QUICK_GPU_TARGET_SUFFIX}.MPI $<TARGET_OBJECTS:sander_base_obj_mpi> ${RISM_SANDER_IFACE_MPI})
            target_link_libraries(sander.quick.${QUICK_GPU_TARGET_SUFFIX}.MPI mpi_fortran)
            list(APPEND CUDA_MPI_EXECUTABLES sander.quick.${QUICK_GPU_TARGET_SUFFIX}.MPI)
        endif()

    endif()

endif()

# Things that apply to all executables and/or libraries
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# force linker language to Fortran in case of auto-detection failure
set_property(TARGET ${EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${OPENMP_EXECUTABLES} ${LIBRARIES} PROPERTY LINKER_LANGUAGE Fortran)

# Serial libraries
targets_link_libraries(${EXECUTABLES} ${CUDA_EXECUTABLES} ${OPENMP_EXECUTABLES} ${LIBRARIES} LIBRARIES emil sff_fortran pbsa_sander sebomd netlib amber_common netcdff netcdf ${KMMD_LIB})
if(USE_FFT)
	targets_link_libraries(${EXECUTABLES} ${CUDA_EXECUTABLES} ${OPENMP_EXECUTABLES} ${LIBRARIES} LIBRARIES rism fftw)
endif()

# MPI libraries
if(MPI)
	targets_link_libraries(${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} LIBRARIES emil_mpi sff_fortran_mpi libsqm_mpi pbsa_sander_mpi sebomd_mpi netlib amber_common netcdff netcdf ${KMMD_LIB})

    if(USE_FFT)
        targets_link_libraries(${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} LIBRARIES rism_mpi fftw fftw_mpi)
    endif()
endif()

# QUICK
if(BUILD_QUICK)
    # Link to appropriate QUICK lib
    targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${LIBRARIES} LIBRARIES libquick libcew)
    add_dependencies(sander_base_obj libquick libcew)

    if(MPI)
        targets_link_libraries(${MPI_EXECUTABLES} LIBRARIES libquick_mpi libcew_mpi)
        add_dependencies(sander_base_obj_mpi libquick_mpi libcew_mpi)
    endif()

    if(OPENMP)
        add_dependencies(sander_base_obj_omp libquick libcew)
    endif()

    if(CUDA OR HIP)
        targets_link_libraries(${CUDA_EXECUTABLES} LIBRARIES libquick_${QUICK_GPU_TARGET_SUFFIX} libcew)
        if(MPI)
            targets_link_libraries(${CUDA_MPI_EXECUTABLES} LIBRARIES  libquick_mpi_${QUICK_GPU_TARGET_SUFFIX} libcew_mpi)
        endif()
    endif()
endif()

# TCPB
if(BUILD_TCPB)
    targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${LIBRARIES} LIBRARIES libtcpb)
    add_dependencies(sander_base_obj libtcpb)

    if(MPI)
        targets_link_libraries(${MPI_EXECUTABLES} LIBRARIES libtcpb)
        add_dependencies(sander_base_obj_mpi libtcpb)
    endif()

    if(OPENMP)
        add_dependencies(sander_base_obj_omp libtcpb)
    endif()

    if(CUDA OR HIP AND BUILD_QUICK)
        targets_link_libraries(${CUDA_EXECUTABLES} LIBRARIES libtcpb)
        if(MPI)
            targets_link_libraries(${CUDA_MPI_EXECUTABLES} LIBRARIES libtcpb)
        endif()
    endif()
endif()

# reaxff-puremd
if(BUILD_REAXFF_PUREMD)
    # ${EXECUTABLES} always contains at least the sander executable
    targets_link_libraries(${EXECUTABLES} LIBRARIES libspuremd)

    if(NOT LIBRARIES STREQUAL "")
        targets_link_libraries(${LIBRARIES} LIBRARIES libspuremd)
    endif()

    if(OPENMP AND NOT OPENMP_EXECUTABLES STREQUAL "")
        targets_link_libraries(${OPENMP_EXECUTABLES} LIBRARIES libspuremd)
    endif()

    if(MPI)
        targets_link_libraries(${MPI_EXECUTABLES} LIBRARIES libspuremd)
    endif()

    if(CUDA AND NOT CUDA_EXECUTABLES STREQUAL "")
        targets_link_libraries(${CUDA_EXECUTABLES} LIBRARIES libspuremd)
	if(MPI AND NOT CUDA_MPI_EXECUTABLES STREQUAL "")
            targets_link_libraries(${CUDA_MPI_EXECUTABLES} LIBRARIES libspuremd)
        endif()
    endif()
endif()


if(dftbplus_EXTERNAL)
  #find_package(DftbPlus REQUIRED)
  get_target_property(DFTBPLUS_INC_DIR DftbPlus::dftbplus INTERFACE_INCLUDE_DIRECTORIES)
  include_directories(${DFTBPLUS_INC_DIR})
  targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES DftbPlus::DftbPlus)
  add_compile_definitions(DFTBPLUS)
  
  get_target_property(DFTBPLUS_LINK_LIBRARIES DftbPlus::dftbplus INTERFACE_LINK_LIBRARIES)
  if("s-dftd3::s-dftd3" IN_LIST DFTBPLUS_LINK_LIBRARIES)
    add_compile_definitions(WITH_SDFTD3)
    #if(xtb_EXTERNAL)
    #  message(FATAL_ERROR "Cannot link xtb and dftbplus into sander at the same time because dftbplus was configured with the -DWITH_SDFT3=TRUE cmake option. Either you need to reconfigure and reinstall dftbplus without said option, or you need to reconfigure AmberTools with -DUSE_DFTBPLUS=FALSE or -DUSE_XTB=FALSE")
    #endif()
  endif()
endif()



if(xtb_EXTERNAL)
    # find_package(XTB REQUIRED)
    # ${EXECUTABLES} always contains at least the sander executable
    targets_link_libraries(${EXECUTABLES} sander_base_obj LIBRARIES xtb::xtb)
    
    add_compile_definitions(XTB)

    if(NOT LIBRARIES STREQUAL "")
        targets_link_libraries(${LIBRARIES} LIBRARIES xtb::xtb)
    endif()

    if(OPENMP AND NOT OPENMP_EXECUTABLES STREQUAL "")
        targets_link_libraries(${OPENMP_EXECUTABLES} sander_base_obj_omp LIBRARIES xtb::xtb)
    endif()

    if(MPI)
        targets_link_libraries(${MPI_EXECUTABLES} sander_base_obj_mpi LIBRARIES xtb::xtb)
    endif()

    if(CUDA AND NOT CUDA_EXECUTABLES STREQUAL "")
        targets_link_libraries(${CUDA_EXECUTABLES} LIBRARIES xtb::xtb)
	    if(MPI AND NOT CUDA_MPI_EXECUTABLES STREQUAL "")
            targets_link_libraries(${CUDA_MPI_EXECUTABLES} LIBRARIES xtb::xtb)
        endif()
    endif()
endif()


if(deepmd-kit_EXTERNAL)
  get_target_property(DeePMD_INC_DIR DeePMD::deepmd_c INTERFACE_INCLUDE_DIRECTORIES)
  include_directories(${DeePMD_INC_DIR})
  targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES DeePMD::deepmd_c)
  add_compile_definitions(DPRC)
endif()
  

# Boost headers
targets_link_libraries(${EXECUTABLES} ${LIBRARIES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} LIBRARIES boost_headers)

targets_link_libraries(${EXECUTABLES} ${LIBRARIES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} LIBRARIES libsqm)

if(PLUMED_RUNTIME_LINK)
	targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES dl)
	set_target_properties(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} PROPERTIES ENABLE_EXPORTS TRUE)
else()
	if(plumed_ENABLED)
		targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES plumed::plumed)
	endif()
endif()

if(mbx_ENABLED)
    targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES MBX::mbx)
endif()

if(lio_ENABLED)
    targets_link_libraries(${EXECUTABLES} ${OPENMP_EXECUTABLES} ${MPI_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} ${CUDA_EXECUTABLES} ${LIBRARIES} LIBRARIES lio::lio)
endif()

install(TARGETS ${EXECUTABLES} DESTINATION ${BINDIR})
install(TARGETS ${MPI_EXECUTABLES} DESTINATION ${BINDIR} COMPONENT MPI)
install(TARGETS ${CUDA_EXECUTABLES} ${CUDA_MPI_EXECUTABLES} DESTINATION ${BINDIR} COMPONENT CUDA)
install(TARGETS ${OPENMP_EXECUTABLES} DESTINATION ${BINDIR} COMPONENT OpenMP)
install_libraries(${LIBRARIES})
