
file(GLOB SOURCES *.cpp *.hpp)
add_executable(ndfes-path ${SOURCES})
target_link_libraries(ndfes-path PRIVATE xmlio ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES})
set_target_properties(ndfes-path PROPERTIES 
  OUTPUT_NAME "ndfes-path"
  INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
)

install(TARGETS ndfes-path DESTINATION ${CMAKE_INSTALL_BINDIR})

if (OpenMP_FOUND)
  add_executable(ndfes-path${OMP_SUFFIX} ${SOURCES})
  target_link_libraries(ndfes-path${OMP_SUFFIX} PRIVATE xmlio ${BLAS_LIBRARIES} ${NLOPT_LIBRARIES} OpenMP::OpenMP_CXX)
  set_target_properties(ndfes-path${OMP_SUFFIX} PROPERTIES 
    OUTPUT_NAME "ndfes-path${OMP_SUFFIX}"
    INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
  )
  target_compile_definitions(ndfes-path${OMP_SUFFIX} PRIVATE WITH_OPENMP)
  install(TARGETS ndfes-path${OMP_SUFFIX} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
