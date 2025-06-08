
file(GLOB SOURCES *.cpp *.hpp)
add_executable(edgembarcli ${SOURCES})
target_link_libraries(edgembarcli PRIVATE edgembarlib)
set_target_properties(edgembarcli PROPERTIES 
  OUTPUT_NAME "edgembar"
  INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
)

install(TARGETS edgembarcli DESTINATION ${CMAKE_INSTALL_BINDIR})

if (OpenMP_FOUND)
  add_executable(edgembarcli${OMP_SUFFIX} ${SOURCES})
  target_link_libraries(edgembarcli${OMP_SUFFIX} PRIVATE edgembarlib${OMP_SUFFIX} OpenMP::OpenMP_CXX)
  set_target_properties(edgembarcli${OMP_SUFFIX} PROPERTIES 
    OUTPUT_NAME "edgembar${OMP_SUFFIX}"
    INSTALL_RPATH \$ORIGIN/../lib:\$ORIGIN/../lib64
    )
  target_compile_definitions(edgembarcli${OMP_SUFFIX} PRIVATE WITH_OPENMP)
  install(TARGETS edgembarcli${OMP_SUFFIX} DESTINATION ${CMAKE_INSTALL_BINDIR})
endif()
