project(xtalutil C)


#Sources for the libXtal library
set(XTAL_SOURCES crdmanip.c grid.c matrix.c
    pdbRead.c rand2.c topRead.c vector.c)

#------------------------------------------------------------------------------------------

add_compile_options(${OPT_CFLAGS})

add_library(Xtal ${XTAL_SOURCES})
target_link_libraries(Xtal C::Math)
install_libraries(Xtal)

#------------------------------------------------------------------------------------------

#since each program is a single file with the same dependencies, we list them here and use a loop to create them

set(XTAL_PROGRAMS AddToBox ChBox PropPDB UnitCell)

#------------------------------------------------------------------------------------------

#since libxtal is built optimized, I assume that these should be optimized as well
#the original makefiledidn't use any optimization flags.
add_compile_options(${OPT_CFLAGS})

foreach(XTAL_PROG ${XTAL_PROGRAMS})
    add_executable(${XTAL_PROG} "${XTAL_PROG}.c")
    target_link_libraries(${XTAL_PROG} Xtal)
endforeach()

install(TARGETS ${XTAL_PROGRAMS} DESTINATION ${BINDIR})
