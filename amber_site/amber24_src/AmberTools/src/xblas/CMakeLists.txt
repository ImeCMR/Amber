#AMBER buildfile for xblas
project(xblas)

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(XBLAS_LIBRARY ${BUILD_DIR}/libxblas-amb.a)

list_to_space_separated(CMAKE_EXE_LINKER_FLAGS_SPC ${CMAKE_EXE_LINKER_FLAGS})

# Determine C flags
# --------------------------------------------------------------------

set(XBLAS_CFLAGS "${CMAKE_C_FLAGS} ${OPT_CFLAGS_SPC}")

if(SHARED)
	#we're kind of breaking the rules here by using this internal CMake variable, 
	#but we need to pass the position independent code flag
	list_to_space_separated(SHARED_FLAGS ${CMAKE_SHARED_LIBRARY_C_FLAGS})
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} ${SHARED_FLAGS}")
endif()

# Fortran name mangling
include(FortranCInterface)
FortranCInterface_VERIFY() #check that Fortran can link to C

if("${FortranCInterface_GLOBAL_CASE}" STREQUAL UPPER)
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} -DCONFIG_FC_UCASE")
endif()

if("${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL _)
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} -DCONFIG_FC_UNDERSCORE")
elseif("${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL __)
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} -DCONFIG_FC_DBL_UNDERSCORE")
endif()

# check whether FPU control code is needed
check_all_includes(C fpu_control.h)

if(HAVE_FPU_CONTROL_H)
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} -DHAVE_FPU_CONTROL_H")
endif()

if("${TARGET_ARCH}" STREQUAL "x86_64" OR "${TARGET_ARCH}" STREQUAL "ia64" OR "${TARGET_ARCH}" STREQUAL "i386")
	set(XBLAS_CFLAGS "${XBLAS_CFLAGS} -Dx86")
endif()

# Link Libraries
# --------------------------------------------------------------------

set(XBLAS_EXTRA_LIBS ${CMath_LIBRARIES}) # used in configuration scripts

# M4 configuration
# --------------------------------------------------------------------

set(XBLAS_M4_OPTS "")

# Non-GNU m4 require expanded buffer size
# Detect GNU m4 by its behaviour on __gnu__.
# gnu m4 will output "__gnu__" with an empty string, and non-gnu will output it verbatim.

file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/m4check.in "__gnu__
")

execute_process(COMMAND ${M4}
	RESULT_VARIABLE M4_RESULT
	OUTPUT_VARIABLE M4_OUTPUT
	ERROR_VARIABLE M4_ERROR
	INPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/m4check.in)

if(NOT "${M4_RESULT}" EQUAL 0)
	message(FATAL_ERROR "Failed to check vendor of m4.  ${M4} failed with exit code ${M4_RESULT} and output ${M4_ERROR}")
endif()

if("${M4_OUTPUT}" MATCHES "__gnu__")
	# non-GNU M4
	set(XBLAS_M4_OPTS "${XBLAS_M4_OPTS} -B65536")
endif()

# hardcoded config options
#build plain (no extended or mixed precision) BLAS
set(ENABLE_PLAIN_BLAS FALSE)
set(ENABLE_FORTRAN FALSE)

if(NOT ENABLE_PLAIN_BLAS)
	set(XBLAS_M4_OPTS "${XBLAS_M4_OPTS} -D no_plain_blas")
endif()

if(NOT ENABLE_FORTRAN)
	set(XBLAS_M4_OPTS "${XBLAS_M4_OPTS} -D no_f2c")
endif()

# --------------------------------------------------------------------

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/make.inc.in.cmake ${CMAKE_CURRENT_BINARY_DIR}/make.inc)

# --------------------------------------------------------------------

#note: XBLAS does not support out-of-source builds, so we copy the entire source directory to the build tree

ExternalProject_Add(xblas_build
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}
	SOURCE_DIR ${BUILD_DIR}
    BUILD_IN_SOURCE TRUE
	DOWNLOAD_DIR ${BUILD_DIR}
	DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${BUILD_DIR} 
	COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/make.inc ${BUILD_DIR}/make.inc
	CONFIGURE_COMMAND ""
	BUILD_COMMAND ${MAKE_COMMAND} lib-amb
	INSTALL_COMMAND ""
	BUILD_BYPRODUCTS ${XBLAS_LIBRARY}) #no install command

add_library(xblas STATIC IMPORTED GLOBAL)
add_dependencies(xblas xblas_build)
set_property(TARGET xblas PROPERTY IMPORTED_LOCATION ${XBLAS_LIBRARY})

install(FILES ${XBLAS_LIBRARY} DESTINATION ${LIBDIR})

