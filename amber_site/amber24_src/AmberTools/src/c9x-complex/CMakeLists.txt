#AMBER CMake buildfile for c9x-complex
project(c9x-complex C)

include(TestBigEndian)
test_big_endian(BIGENDIAN_BOOL)
#the header expects 1 or 0 for the BIGENDIAN constant
if(${BIGENDIAN_BOOL})
	set(BIGENDIAN_VALUE 1)
else()
	set(BIGENDIAN_VALUE 0)
endif() 

#set up the configure file
set(CONFIG_HEADER_PATH "${CMAKE_CURRENT_BINARY_DIR}/mconf.h")  
configure_file(mconf.h.in ${CONFIG_HEADER_PATH})

#append the do-not-edit warning
file(READ ${CONFIG_HEADER_PATH} CONFIG_HEADER_CONTENTS)
set(CONFIG_HEADER_CONTENTS "/*-------------------------------------------------------------------------------------------------------------------------------------------------
 * DO NOT EDIT!!!! This file is generated from mconf.h.in in the source tree.  Change these defines there, not here!
 -------------------------------------------------------------------------------------------------------------------------------------------------*/

${CONFIG_HEADER_CONTENTS}")
file(WRITE ${CONFIG_HEADER_PATH} "${CONFIG_HEADER_CONTENTS}")

#setup the libmc library
add_library(mc STATIC cmplx.c clog.c cgamma.c stubs.c)

target_include_directories(mc PRIVATE ${CONFIG_HEADER_DIR} PUBLIC .)

