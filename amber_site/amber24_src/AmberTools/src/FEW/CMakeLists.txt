#CMake buildfile for FEW
project(FEW)

# draw_membrane2 executable
add_compile_options(${NO_OPT_CFLAGS})
add_executable(draw_membrane2 miscellaneous/draw_membrane2.c)

target_link_libraries(draw_membrane2 C::Math)

install(TARGETS draw_membrane2 DESTINATION ${BINDIR})

install(PROGRAMS FEW.pl DESTINATION ${BINDIR})


#install the dependency perl modules
add_subdirectory(additional_libs)
add_subdirectory(libs)
add_subdirectory(examples)