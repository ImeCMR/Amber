project(xray Fortran)

#list of all programs built by this file
#used for setting optimization flags
set(F_PROGRAMS add_pdb add_xray)

#shell scripts to install
set(SHELL_PROGS XrayPrep)
#------------------------------------------------------------------------------

add_executable(add_pdb add_pdb.F90)
target_link_libraries(add_pdb amber_common)

add_executable(add_xray add_xray.F90)
target_link_libraries(add_xray amber_common)

#optimization flags
set_property(TARGET ${F_PROGRAMS} PROPERTY COMPILE_OPTIONS ${NO_OPT_FFLAGS})

#install everything
#------------------------------------------------------------------------------

install(TARGETS ${F_PROGRAMS} DESTINATION ${BINDIR})
install(PROGRAMS ${SHELL_PROGS} DESTINATION ${BINDIR})

