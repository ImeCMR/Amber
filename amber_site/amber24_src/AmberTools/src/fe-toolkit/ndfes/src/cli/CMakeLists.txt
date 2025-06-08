
file(GLOB SOURCES *.cpp *.hpp)
add_executable(ndfes ${SOURCES})
target_link_libraries(ndfes PRIVATE ndfeslib)
set_target_properties(ndfes PROPERTIES 
  OUTPUT_NAME "ndfes"
  INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
)

install(TARGETS ndfes DESTINATION ${CMAKE_INSTALL_BINDIR})

if (OpenMP_FOUND)
  add_executable(ndfes${OMP_SUFFIX} ${SOURCES})
  target_link_libraries(ndfes${OMP_SUFFIX} PRIVATE ndfeslib${OMP_SUFFIX} OpenMP::OpenMP_CXX)
  set_target_properties(ndfes${OMP_SUFFIX} PROPERTIES 
    OUTPUT_NAME "ndfes${OMP_SUFFIX}"
    INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
  )
  target_compile_definitions(ndfes${OMP_SUFFIX} PRIVATE WITH_OPENMP)
  install(TARGETS ndfes${OMP_SUFFIX} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
