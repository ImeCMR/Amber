# check if we are building as part of AmberTools
# PROJECT_NAME will be overriden by `project` command,
# so we need to check it before calling `project`
if (NOT DEFINED IS_AMBERTOOLS)
  if (PROJECT_NAME STREQUAL "Amber")
    set(IS_AMBERTOOLS TRUE)
  else()
    set(IS_AMBERTOOLS FALSE)
  endif()
endif()

#
# Let's use 3.10, if possible, because that's the highest
# requirement of any part of AmberTools
#
if (IS_AMBERTOOLS)
  cmake_minimum_required(VERSION 3.10)
else()
  cmake_minimum_required(VERSION 3.12)
endif()

set(FETOOLKITVERSION 3.1)

#
# Fortran is needed so we get -lgfortran when linking to some
# blas/lapack implementations
#
project(fe-toolkit VERSION "${FETOOLKITVERSION}" LANGUAGES CXX Fortran)


add_subdirectory(edgembar)
add_subdirectory(ndfes)
  

