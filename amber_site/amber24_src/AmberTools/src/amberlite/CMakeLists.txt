project(amberlite C)

#I assume that these programs are supposed to be optimized
add_compile_options(${OPT_CFLAGS})

#add the three nab executables
add_nab_executable(ffgbsa NAB_SOURCES ffgbsa.nab)
add_nab_executable(minab NAB_SOURCES minab.nab)
add_nab_executable(mdnab NAB_SOURCES mdnab.nab)

#install them
install(TARGETS ffgbsa minab mdnab DESTINATION ${BINDIR})

#install the python library
install_python_library()