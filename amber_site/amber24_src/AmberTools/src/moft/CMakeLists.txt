set(METATWIST_SOURCES metatwist.cpp)

add_executable(metatwist ${METATWIST_SOURCES})

target_include_directories(metatwist PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(metatwist PRIVATE ${OPT_CXXFLAGS})
target_link_libraries(metatwist fftw boost_program_options boost_iostreams boost_regex boost_system boost_timer boost_chrono boost_filesystem boost_graph)

install(TARGETS metatwist DESTINATION ${BINDIR})