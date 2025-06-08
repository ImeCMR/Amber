#tree CMake file for AMBER

foreach(TOOL ${AMBER_TOOLS})
	if("${TOOL}" IN_LIST TOOLS_IN_SRC)
		add_subdirectory(${TOOL})
	endif()
endforeach()
