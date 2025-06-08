project(pnetcdf C Fortran)

# Hardcoded settings
# --------------------------------------------------------------------

set(ENABLE_FORTRAN TRUE)

if(NOT MPI)
	message(FATAL_ERROR "This directory should not be built with MPI turned off!")
endif()

set(PACKAGE_NAME pnetcdf)
set(PACKAGE_VERSION ?)
set(PACKAGE_BUGREPORT "")
set(PACKAGE_STRING ${PACKAGE_NAME}-${PACKAGE_VERSION})

# Fortran and C linking test
# --------------------------------------------------------------------
include(FortranCInterface)
FortranCInterface_VERIFY() #check that Fortran can link to C

test(F77_NAME_LOWER "${FortranCInterface_GLOBAL_CASE}" STREQUAL LOWER)
test(F77_NAME_UPPER "${FortranCInterface_GLOBAL_CASE}" STREQUAL UPPER)
test(F77_NAME_LOWER_2USCORE F77_NAME_LOWER AND ("${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL __))
test(F77_NAME_LOWER_USCORE F77_NAME_LOWER AND ("${FortranCInterface_GLOBAL_SUFFIX}" STREQUAL _))

# Set up MPI for configure tests and code compilation
# --------------------------------------------------------------------
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})

set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})

# Function checks
# --------------------------------------------------------------------
check_symbol_exists(alloca alloca.h HAVE_ALLOCA)

check_symbol_exists(MPI_MODE_RDWR mpi.h HAVE_MPI_IO)
if(NOT HAVE_MPI_IO)
	message(FATAL_ERROR "Your mpi implementation does not support MPI-IO! Internal pnetcdf cannot be built!")
endif()

# Check for presence of various MPI_COMBINERs.  These could be enums, so we have to do symbol checks.
foreach(COMBINER_NAME DUP HVECTOR_INTEGER HINDEXED_INTEGER SUBARRAY DARRAY RESIZED STRUCT_INTEGER INDEXED_BLOCK F90_REAL F90_INTEGER F90_COMPLEX)
	check_all_symbols(mpi.h "MPI_COMBINER_${COMBINER_NAME}")
endforeach()

# Check for presence of Fortran types
#
# We do this for a couple of reasons.  First, the MPI might have been
# built without Fortran support, in which case these types might not
# exist.  Second, we need to map these types to corresponding C types
# where possible to simplify processing at runtime.
check_all_symbols(mpi.h 
MPI_CHARACTER
MPI_REAL
MPI_INTEGER
MPI_DOUBLE_PRECISION
MPI_INTEGER1
MPI_INTEGER2
MPI_INTEGER4
MPI_INTEGER8
MPI_INTEGER16
MPI_REAL4
MPI_REAL8
MPI_REAL16
MPI_COMPLEX8
MPI_COMPLEX16
MPI_COMPLEX32
MPI_UB
MPI_LB)

check_all_functions(strerror MPI_Request_get_status MPI_Type_dup MPI_Info_dup)

test(NO_STRERROR NOT HAVE_STRERROR)

# the nonblocking routines need these functions
test(ENABLE_NONBLOCKING HAVE_MPI_REQUEST_GET_STATUS AND HAVE_MPI_TYPE_DUP)

# Type size checks
# --------------------------------------------------------------------
check_all_types(ptrdiff_t ssize_t uchar MPI_Offset short int long float double size_t off_t "((struct stat*)0)->st_blksize")

include(${CMAKE_CURRENT_SOURCE_DIR}/../netcdf-fortran-4.4.4/MatchNetcdfFortranTypes.cmake)

# I feel like we can assume float is IEEE compliant
test(NO_IEEE_FLOAT NOT HAVE_FLOAT)

# Header checks
# --------------------------------------------------------------------
check_all_includes(C alloca.h stdlib.h sys/types.h mpio.h)

# Configure header
# --------------------------------------------------------------------
configure_file(ncconfig.cmake.in.h ${CMAKE_CURRENT_BINARY_DIR}/ncconfig.h)
configure_file(nfconfig.cmake.in.inc ${CMAKE_CURRENT_BINARY_DIR}/nfconfig.inc)

include_directories(${CMAKE_CURRENT_BINARY_DIR})

add_subdirectory(src)

if(ENABLE_FORTRAN)
	set(PNETCDF_FORTRAN_OBJ $<TARGET_OBJECTS:pnetcdf_fortran_obj>)
else()
	set(PNETCDF_FORTRAN_OBJ "")
endif()

add_library(pnetcdf $<TARGET_OBJECTS:pnetcdf_c_obj> ${PNETCDF_FORTRAN_OBJ})

target_link_libraries(pnetcdf mpi_c C::Math)

target_include_directories(pnetcdf INTERFACE src/lib)
install_libraries(pnetcdf)