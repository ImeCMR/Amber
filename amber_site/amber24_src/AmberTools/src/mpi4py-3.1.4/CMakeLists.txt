# CMake buildfile for mpi4py
# This is kind of complicated because MPI has all sorts of flags and paths that need to be set.
# mpi4py handles this using a configuration file, mpi.cfg, which has entries for all of the different types of
# MPI.  CMake's FindMPI module gets us this information easily, so all we have to do is proxy it to the configuration file.
# We do NOT just pass the MPI wrappers to mpi4py, because not all MPI implementations have wrappers

project(mpi4py)

set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/mpi4py-build)

set(STAMP_FILE ${BUILD_DIR}/mpi4py-build.stamp)

set_property(DIRECTORY PROPERTY ADDITIONAL_MAKE_CLEAN_FILES ${BUILD_DIR})

#find all python scripts and C files as dependencies
file(GLOB_RECURSE MPI4PY_SOURCES "*.py" "*.c" "*.pxd" "*.h")

# --------------------------------------------------------------------
# Get the MPI flags into the config file

# separator for lists in mpi.cfg
if(HOST_WINDOWS)
	set(LIST_SEP ";")
else()
	set(LIST_SEP ":")
endif()

string(REPLACE ";" "${LIST_SEP}" MPI4PY_INCDIRS "${MPI_C_INCLUDE_PATH}")

# MPI headers might also live in the prefix path (as with MS-MPI), in which case mpi4py won't be able to find them
if(DEFINED CMAKE_PREFIX_PATH)
	if(NOT "${CMAKE_PREFIX_PATH}" STREQUAL "")
		set(MPI4PY_INCDIRS "${MPI4PY_INCDIRS}${LIST_SEP}${CMAKE_PREFIX_PATH}/include")
	endif()
endif()

set(MPI4PY_LIBS_LIST "")
set(MPI4PY_LIBDIR_LIST "")

# split the link libraries into names and paths
foreach(LIB ${MPI_C_LIBRARIES})
	get_filename_component(LIB_FOLDER ${LIB} PATH)
	get_filename_component(LIB_FILENAME ${LIB} NAME_WE)
	
	string(REPLACE "lib" "" LIB_NAME ${LIB_FILENAME})
	
	list(APPEND MPI4PY_LIBS_LIST ${LIB_NAME})
	list(APPEND MPI4PY_LIBDIR_LIST ${LIB_FOLDER})
endforeach()

if(TARGET_OSX)
	set(MPI4PY_LINK_FLAGS "${MPI_C_LINK_FLAGS} -Wl,-rpath,@loader_path/../../..")
else()
	set(MPI4PY_LINK_FLAGS ${MPI_C_LINK_FLAGS})
endif()

# most of the libraries will probably be in the link folder
list(REMOVE_DUPLICATES MPI4PY_LIBDIR_LIST)

list_to_space_separated(MPI4PY_LIBS ${MPI4PY_LIBS_LIST})

string(REPLACE ";" "${LIST_SEP}" MPI4PY_LIBDIRS "${MPI4PY_LIBDIR_LIST}")

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/mpi.cfg.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/mpi.cfg)

# --------------------------------------------------------------------
# build command

add_custom_command(OUTPUT ${STAMP_FILE}
    COMMAND CC=${MPI_C_COMPILER} ${PYTHON_EXECUTABLE} setup.py build ${PYTHON_COMPILER_ARG} -b ${BUILD_DIR} "--mpi=cmake-mpi,${CMAKE_CURRENT_BINARY_DIR}/mpi.cfg"
	COMMAND ${CMAKE_COMMAND} -E touch ${STAMP_FILE}
	WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
	DEPENDS ${MPI4PY_SOURCES} 
	VERBATIM
	COMMENT "Building mpi4py")


#We want to build the python library during the build step so as to catch any build errors
add_custom_target(mpi4py ALL DEPENDS ${STAMP_FILE})

install(CODE "
	${FIX_BACKSLASHES_CMD}
	execute_process(
        COMMAND  \"${PYTHON_EXECUTABLE}\" 
    	\"${CMAKE_CURRENT_SOURCE_DIR}/setup.py\" 
    	 build -b \"${BUILD_DIR}\" ${PYTHON_COMPILER_ARG} \"--mpi=cmake-mpi,${CMAKE_CURRENT_BINARY_DIR}/mpi.cfg\"
    	 install -f ${PYTHON_PREFIX_ARG} --single-version-externally-managed --root /
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})" COMPONENT Python)
