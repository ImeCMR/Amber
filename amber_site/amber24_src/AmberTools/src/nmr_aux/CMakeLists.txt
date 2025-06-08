#AMBER buildfile for nmr_aux
project(nmr_aux Fortran C)

# add_subdirectory(fantasian)
add_subdirectory(prepare_input)

#this directory seemed to have been disabled in the makefile, not sure if it should be built
#add_subdirectory(rotdif)

install(PROGRAMS senergy sviol sviol2 DESTINATION ${BINDIR})
