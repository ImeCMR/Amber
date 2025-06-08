project(reduce CXX)

# platform checks

# Require C++-11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_STANDARD_REQUIRED ON)
check_all_functions(sscanf)

if(HAVE_SSCANF)
	add_definitions(-DHASSSCANFEXTERN)
endif()

# I think we can assume this at this point...
add_definitions(-DBOOLPREDEFINED)

# we tried template exporting, and it didn't work.
add_definitions(-DINCTEMPLATEDEFNS)

# template operators are well supported
add_definitions(-DBRACKETOPERPARMS)
# --------------------------------------------------------------------

add_subdirectory(libpdb)
add_subdirectory(toolclasses)
add_subdirectory(reduce_src)

install(FILES reduce_wwPDB_het_dict.txt DESTINATION ${DATADIR})
