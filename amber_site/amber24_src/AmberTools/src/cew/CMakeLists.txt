cmake_minimum_required(VERSION 3.8.0)
project(cew NONE)

#version number
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
set(cew_MAJOR_VERSION 1)
set(cew_MINOR_VERSION 0)
set(cew_TWEAK_VERSION 0)

set(cew_VERSION "${cew_MAJOR_VERSION}.${cew_MINOR_VERSION}.${cew_TWEAK_VERSION}")

if(NOT INSIDE_AMBER)
	# initialization and include paths
	#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
	include(cmake/AmberBuildSystemInit.cmake)

	#now enable the languages
	enable_language(C CXX Fortran)

	include(AmberBuildSystem2ndInit)

	set(LINALG_LIBS_REQUIRED FALSE)
	set(NEEDED_3RDPARTY_TOOLS )
	set(REQUIRED_3RDPARTY_TOOLS fftw)
	include(3rdPartyTools)

	include(CompilerFlags)

	#CPack setup
	# --------------------------------------------------------------------
	set(PACKAGE_NAME "cew")
	set(PACKAGE_FILENAME "cew")

	set(BUNDLE_IDENTIFIER org.yorklab.cew)
	set(BUNDLE_SIGNATURE cew)
	include(Packaging)

	# build options
	option(INSTALL_HEADERS "Copy headers to the include/cew folder of the install directory." TRUE)

	# Ensure environment is read correctly
	string(REPLACE ":" ";" _lib_path "$ENV{LD_LIBRARY_PATH}")
	set( CMAKE_LIBRARY_PATH ${_lib_path} )

	string(REPLACE ":" ";" _inc_path "$ENV{CPATH}")
	set( CMAKE_INCLUDE_PATH ${_inc_path} )
else()
	set(INSTALL_HEADERS TRUE)
endif()

# on Windows, make MSVC auto-create import libraries just like MinGW does
set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS TRUE)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Parse cewfiles (this file contains a list of all files to be compiled and the respective objects)

# read each non-empty line into an element of a list
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/cewfiles CEWFILES_CONTENTS)

# get rid of backslashes
string(REPLACE "\\" "" CEWFILES_CONTENTS "${CEWFILES_CONTENTS}")

# name of list that we are currently appending to
set(LIST_NAME "")

foreach(LINE ${CEWFILES_CONTENTS})

	# ignore comment lines
	if(NOT "${LINE}" MATCHES "^#")

		# extract the name of the source file mentioned in the line (a string after whitespace or an equals sign)
		string(REGEX MATCH "[^ :=]+\.(o|cpp|c|LIBCEW\.o)" SOURCE_FILE_NAME "${LINE}")

		# get name of variable that the following list is being set to
		# must exclude parentheses so that we don't match dereferences of other variables
		string(REGEX MATCH "[^$\(\)]+=" VARIABLE_NAME "${LINE}")

		# if we are starting a new source list, update LIST_NAME accordingly
		if(NOT "${VARIABLE_NAME}" STREQUAL "")
			string(REPLACE "=" "" VARIABLE_NAME "${VARIABLE_NAME}")
			set(LIST_NAME ${VARIABLE_NAME})
		endif()

		# did we get a new source file?
		if(NOT "${SOURCE_FILE_NAME}" STREQUAL "")

			if("${LIST_NAME}" STREQUAL "")
				message(FATAL_ERROR "cewfiles parser error: got source files before any source lists!")
			endif()

			# get rid of LIBCEW.o suffix if it exists
			string(REPLACE "LIBCEW.o" "cpp" SOURCE_FILE_NAME "${SOURCE_FILE_NAME}")

			# also convert .o to .cpp (used in some variables)
			string(REPLACE ".o" ".cpp" SOURCE_FILE_NAME "${SOURCE_FILE_NAME}")


			list(APPEND ${LIST_NAME} ${SOURCE_FILE_NAME})
		endif()

		#message("\"${LINE}\" - SFN: \"${SOURCE_FILE_NAME}\" - VN: \"${VARIABLE_NAME}\"")
	endif()
endforeach()

# The above loop will create the following variables:
# SOURCES - C++ sources for libcew
# OBJECTS - C++ sources for libcew

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

# with icc, cew needs -fp-model source in order to produce floating point results that match gcc
set(FP_SOURCE_FLAG "")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
	check_cxx_compiler_flag("-fp-model source" FP_MODEL_SOURCE_WORKS)

	if(FP_MODEL_SOURCE_WORKS)
		set(FP_SOURCE_FLAG "-fp-model source")
	endif()
endif()

# all sources should use optimized compile flags
set_property(SOURCE ${SOURCES} PROPERTY COMPILE_FLAGS "${OPT_CXXFLAGS_SPC} ${OPT_CFLAGS_SPC} ${FP_SOURCE_FLAG}")

#------------------------------------------------------------------------------------------
# libcew library

add_library(libcew ${SOURCES})
target_link_libraries(libcew PUBLIC fftw)
remove_prefix(libcew)
install_libraries(libcew)

if(MPI)
	make_mpi_version(libcew libcew_mpi LANGUAGES CXX)
	target_compile_definitions(libcew_mpi PRIVATE WITH_MPI)
    target_link_libraries(libcew_mpi PUBLIC fftw)
    remove_prefix(libcew_mpi)
    install_libraries(libcew_mpi)
endif()

# Create the headers inside the build folder, just to serve to other tools
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/cew)
file(GLOB CEW_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
add_custom_command(TARGET libcew POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy ${CEW_HEADERS} ${CMAKE_CURRENT_BINARY_DIR}/include/cew)
target_include_directories(libcew INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include/cew)
if(MPI)
    target_include_directories(libcew_mpi INTERFACE ${CMAKE_CURRENT_BINARY_DIR}/include/cew)
endif()

#------------------------------------------------------------------------------------------
# DLL exports/imports

if(SHARED)
	# now set up the definition for other people to use
	target_compile_definitions(libcew INTERFACE CEW_USE_DLL)
    if(MPI)
	    target_compile_definitions(libcew_mpi INTERFACE CEW_USE_DLL)
	endif()
endif()

#------------------------------------------------------------------------------------------
# Header installation

if(INSTALL_HEADERS)

    if(INSIDE_AMBER)
        install(FILES ${CEW_HEADERS} DESTINATION ${AMBERTOOLS_INC_DIR}/cew)
    else()
        install(FILES ${CEW_HEADERS} DESTINATION ${INCDIR}/cew)
    endif()

endif()

#--------------------------------------------------------------
if(NOT INSIDE_AMBER)
	print_build_report()
endif()
