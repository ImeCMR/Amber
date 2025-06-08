#CMake buildfile for ndiff.
#This program is used for checking the output of unit tests.
project(ndiff C)

set(NDIFF_SOURCES ndiff.c store.c awklib.c regexp.c regerror.c)

# --------------------------------------------------------------------

check_all_types("long double")
if(HAVE_LONG_DOUBLE)
	add_definitions(-DHAVE_LONG_DOUBLE)
endif()

# --------------------------------------------------------------------

add_compile_options(${OPT_CFLAGS})

add_executable(ndiff ${NDIFF_SOURCES})