project(sebomd Fortran)
	
set(SEBOMD_SOURCE se_atmchg.F 
         se_backinbox.F90 
         se_balance.F 
         block.F 
         se_bond_order_wrapper.F 
         se_bond_order_write.F90 
         se_bpair.F 
         se_bsort.F 
         se_bspline.F 
         se_clustsub.F 
         se_debug.F 
         se_deeprec.F90
         se_denful.F 
         se_densub.F 
         se_diag.F 
         se_diagp.F 
         se_diat.F 
         se_dihedrpbc.F 
         se_dipole.F 
         se_dmxinter.F 
         se_doferm.F 
         se_doscf.F 
         se_edriver.F 
         se_elecsr.F 
         se_energy.F  
         se_escf.F 
         se_esqr.F 
         se_etimer.F 
         se_fock.F 
         se_fockpme.F 
         se_fourier.F 
         se_fshift.F 
         se_gcart.F 
         se_gdiat.F 
         se_gensub.F 
         se_getpar.F 
         se_glbpnt.F 
         se_goverlp.F 
         se_grepul.F 
         se_ijpair.F 
         se_info_from_sander.F90
         se_initp.F 
         se_mosub.F 
         se_opnfil.F 
         se_overlp.F 
         se_pbcgcres.F 
         se_pbcxyz.F 
         se_pme_calcb.F 
         se_pme_calcchg.F 
         se_pme_calccm1.F 
         se_pme_calcq.F 
         se_pme_calctheta.F 
         se_pme_derec.F 
         se_pme_direct2.F 
         se_pme_qm_direct.F 
         se_pme_qm_recip.F 
         se_pme_recip.F 
         se_pme_setup.F 
         se_pmix.F 
         se_rdall.F 
         se_rdelem.F 
         se_rdnum.F 
         se_rdtail.F 
         se_repul.F 
         se_resclustsub.F 
         se_mpi_allreduce.F 
         se_mpi_vars.F90
         sebomd.F 
         se_corrections_info_from_sander.F90
         se_corrections.F 
         se_init.F90
         se_switch.F90
         se_setbox.F 
         se_setunit.F 
         se_setup.F 
         se_write_charges.F90
         se_sproc.F 
         se_whatis.F 
         se_wrtdmx.F
         se_corrections_params.F90
         se_corrections_module.F90
         se_corrections_tools.F90
         se_inter.F90)

add_compile_options(${OPT_FFLAGS})

if(USE_FFT)
	add_definitions(-DRISMSANDER)
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# must be built as static because it uses sander functions
add_library(sebomd STATIC ${SEBOMD_SOURCE})

config_module_dirs(sebomd ${SEBOMD_MOD_DIR} ${SANDER_COMMON_MOD_DIR})
add_dependencies(sebomd sander_base_obj) #force compilation of sander modules
make_pic_if_needed(sebomd)

target_link_libraries(sebomd netlib amber_common)
install_libraries(sebomd)

if(MPI)
	make_mpi_version(sebomd sebomd_mpi LANGUAGES Fortran INSTALL)
	config_module_dirs(sebomd_mpi ${SEBOMD_MPI_MOD_DIR} ${SANDER_COMMON_MPI_MOD_DIR})
	add_dependencies(sebomd_mpi sander_base_obj_mpi) #force compilation of sander modules
endif()