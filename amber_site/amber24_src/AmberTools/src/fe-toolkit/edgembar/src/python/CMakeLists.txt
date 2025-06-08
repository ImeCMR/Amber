project(edgembar_python)

if (IS_AMBERTOOLS)
    install_python_library()
else()
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --prefix=${CMAKE_INSTALL_PREFIX} ${CMAKE_CURRENT_SOURCE_DIR})")
endif()
