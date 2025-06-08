project(nfe-umbrella-slice CXX)

set(NFE_UMBRELLA_SLICE_SOURCES umbrella.cc 
	dimension.cc utils.cc main.cc)

#------------------------------------------------------------------------------------------

#I'm assuming that this should be optimized
add_compile_options(${OPT_CXXFLAGS})


#------------------------------------------------------------------------------------------

#this is the only header check from the config.h that actually gets used.
check_include_file("getopt.h" HAVE_GETOPT_H)

if(HAVE_GETOPT_H)
	add_definitions(-DHAVE_GETOPT_H)
endif()

#------------------------------------------------------------------------------------------

add_executable(nfe-umbrella-slice ${NFE_UMBRELLA_SLICE_SOURCES})
target_link_libraries(nfe-umbrella-slice netcdf)
install(TARGETS nfe-umbrella-slice DESTINATION ${BINDIR})