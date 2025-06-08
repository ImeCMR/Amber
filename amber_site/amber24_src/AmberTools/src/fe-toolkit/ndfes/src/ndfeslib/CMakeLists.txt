
file(GLOB SOURCES *.cpp *.hpp)
add_library(ndfeslib STATIC ${SOURCES})
target_compile_definitions(ndfeslib PRIVATE PACKAGE_VERSION="${CMAKE_PROJECT_VERSION}")
target_link_libraries(ndfeslib PRIVATE ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES})

if (OpenMP_FOUND)
  add_library(ndfeslib${OMP_SUFFIX} STATIC ${SOURCES})
  target_compile_definitions(ndfeslib${OMP_SUFFIX} PRIVATE PACKAGE_VERSION="${CMAKE_PROJECT_VERSION}")
  target_link_libraries(ndfeslib${OMP_SUFFIX} PRIVATE ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES} OpenMP::OpenMP_CXX)
  target_compile_definitions(ndfeslib${OMP_SUFFIX} PRIVATE WITH_OPENMP)
endif()
