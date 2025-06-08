#CMake ExternalProject makefile for PerlMol

#we can't build in the binary dir because the PerlMol build system generates a makefile which conflicts with CMake's
set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/build)
set(INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/install)

# In some situations, Module::Install seems to think that it was given the INSTALL_BASE argument when it was not
# So, we explicitly blank it
set(CONFIG_COMMAND ${CMAKE_COMMAND} -E env PERL5LIB=${BUILD_DIR}
	${PERL_EXECUTABLE} Makefile.PL 
	PREFIX=${INSTALL_DIR} 
	LIB=${INSTALL_DIR}/${PERL_MODULE_PATH} 
	INSTALL_BASE=
	${PERL_MAKE_SELECT_ARG})


ExternalProject_Add(PerlMol
	PREFIX ${CMAKE_CURRENT_BINARY_DIR}
	SOURCE_DIR ${BUILD_DIR} 
    BUILD_IN_SOURCE TRUE
	DOWNLOAD_DIR ${BUILD_DIR}
	DOWNLOAD_COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${BUILD_DIR}
	CONFIGURE_COMMAND ${CONFIG_COMMAND}
	LIST_SEPARATOR ^
	BUILD_COMMAND ${PERL_MAKE}
	INSTALL_COMMAND  ${PERL_MAKE} install)

install(DIRECTORY ${INSTALL_DIR}/lib DESTINATION ${CMAKE_INSTALL_POSTFIX}.)