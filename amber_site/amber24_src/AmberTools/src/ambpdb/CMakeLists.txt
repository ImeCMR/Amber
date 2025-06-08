set(AMBPDB_SOURCES 
	AmbPDB.cpp 
	DSL.cpp # dummy DataSetList implementation
	Trajout_Single_Stdout.cpp)

add_executable(ambpdb ${AMBPDB_SOURCES})
target_include_directories(ambpdb PUBLIC ../cpptraj/src)
target_link_libraries(ambpdb cpptraj_parm cpptraj_traj cpptraj_file cpptraj_core)
target_compile_options(ambpdb PRIVATE ${OPT_CXXFLAGS})

install(TARGETS ambpdb DESTINATION ${BINDIR})
