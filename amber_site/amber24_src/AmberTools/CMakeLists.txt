#not much to do here
add_subdirectory(src)

if(NOT BUILD_HOST_TOOLS)
	add_subdirectory(test)
endif()

add_subdirectory(examples)
