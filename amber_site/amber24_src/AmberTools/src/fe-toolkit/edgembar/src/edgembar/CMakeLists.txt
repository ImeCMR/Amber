
file(GLOB SOURCES *.cpp *.hpp)
add_library(edgembarlib STATIC ${SOURCES})
target_compile_definitions(edgembarlib PRIVATE PACKAGE_VERSION="${CMAKE_PROJECT_VERSION}")
target_link_libraries(edgembarlib PRIVATE xmlio ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES})

if (OpenMP_FOUND)
  add_library(edgembarlib${OMP_SUFFIX} STATIC ${SOURCES})
  target_compile_definitions(edgembarlib${OMP_SUFFIX} PRIVATE PACKAGE_VERSION="${CMAKE_PROJECT_VERSION}")
  target_link_libraries(edgembarlib${OMP_SUFFIX} PRIVATE xmlio ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES} OpenMP::OpenMP_CXX)
  target_compile_definitions(edgembarlib${OMP_SUFFIX} PRIVATE WITH_OPENMP)  
endif()
