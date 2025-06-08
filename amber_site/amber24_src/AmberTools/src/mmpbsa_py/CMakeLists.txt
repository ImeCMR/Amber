project(mmpbsa_py)


install_python_library()

include_directories(../nabc ../sff)
add_executable(mmpbsa_py_energy mmpbsa_energy.c)
target_link_libraries(mmpbsa_py_energy nabc sff libpbsa)

add_executable(mmpbsa_py_nabnmode mmpbsa_entropy.c)
target_link_libraries(mmpbsa_py_nabnmode nabc sff libpbsa)

install(TARGETS mmpbsa_py_energy mmpbsa_py_nabnmode DESTINATION ${BINDIR} COMPONENT Python)

if(MPI)
	if(NOT mpi4py_ENABLED)
		message(STATUS "WARNING: MMPBSA.py.MPI will not be usable since mpi4py is not enabled")
	endif()
endif()
