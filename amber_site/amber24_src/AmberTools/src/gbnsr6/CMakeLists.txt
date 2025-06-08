#AMBER buildfile for gbnsr6
project(gbnsr6)

set(GBNSR6_SOURCES gbnsr6.F90 pb_read.F90  gb_read.F90 pb_write.F90 getcoor.F90 
	runmd.F90 runmin.F90 force.F90 
	pb_init.F90 sa_driver.F90 np_force.F90 variable_module.F90 
	pb_exmol.F90 NSR6routines.F90  pb_list.F90 timer.F90 
	egb.F90 ene.F90 locmem.F90 myopen.F90 rdparm.F90 
	decomp.F90 rgroup.F90 rfree.F90 debug.F90 
        svdcmp.F90 svbksb.F90 pythag.F90  memory_module.F90 
        gen_dx_file.F90 parms.F90)
        
# --------------------------------------------------------------------

apply_optimization_declarations(${GBNSR6_SOURCES})

add_executable(gbnsr6 ${GBNSR6_SOURCES})
config_module_dirs(gbnsr6 ${GBNSR6_MOD_DIR})

set_property(TARGET gbnsr6 PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(gbnsr6 sff netcdf netlib amber_common)
install(TARGETS gbnsr6 DESTINATION ${BINDIR})


