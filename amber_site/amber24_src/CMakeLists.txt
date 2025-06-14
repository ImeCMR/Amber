#
#           _                      __    ______   ___
#          / \       |\       /|  |  \   |       |   \
#         /   \      | \     / |  |   |  |       |   |
#        /     \     |  \   /  |  |__/   |_____  |__/
#       /-------\    |   \ /   |  |  \   |       |  \
#      /         \   |    v    |  |   |  |       |   \
#     /           \  |         |  |__/   |_____  |    \
#
#
cmake_minimum_required(VERSION 3.8.1)
project(Amber NONE)
#version number
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
set(Amber_MAJOR_VERSION 24)
set(Amber_MINOR_VERSION 0)
set(Amber_PATCH_VERSION 0)

set(Amber_VERSION "${Amber_MAJOR_VERSION}.${Amber_MINOR_VERSION}.${Amber_PATCH_VERSION}")

# flag to submodules to tell them to not run their own version of AmberCmakeConfig
set(INSIDE_AMBER TRUE)


# host tools options
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

option(BUILD_HOST_TOOLS "Build and install only the tools neded to cross-compile Amber.  Disables most build system functionality." FALSE)

option(USE_HOST_TOOLS "Use an installation of Amber Host Tools (pointed to by HOST_TOOLS_DIR) to build amber.  This enables cross compilation." FALSE)

set(HOST_TOOLS_DIR "" CACHE PATH "Path to Amber Host Tools.  Only used if USE_HOST_TOOLS is true.")

if(BUILD_HOST_TOOLS AND USE_HOST_TOOLS)
	message(FATAL_ERROR "Cannot build and use host tools at the same time!  Pick one!")
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Initial configuration

# initialization and include paths
include(cmake/AmberBuildSystemInit.cmake NO_POLICY_SCOPE)

#we don't want to notice src if it is there but empty, or if it's an older version of Amber that lacks CMake support
#so we check for CMakeLists.txt
if(EXISTS "${CMAKE_SOURCE_DIR}/src/CMakeLists.txt")
	message(STATUS "Amber source found, building AmberTools and Amber")
	set(AMBERTOOLS_ONLY FALSE)
else()
	message(STATUS "Amber source not found, only building AmberTools")
	set(AMBERTOOLS_ONLY TRUE)
endif()

# the .gitignore file is stripped out by the release script
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/.gitignore)
	set(AMBER_RELEASE TRUE)
else()
	set(AMBER_RELEASE FALSE)
endif()

if(NOT BUILD_HOST_TOOLS)
	option(INSTALL_TESTS "Whether or not to install ${PROJECT_NAME}'s tests, examples, and benchmarks as part of the main install target.  Be warned, they take up over a gigabyte/" FALSE)
	
	# this option gets passed to install commands
	if(INSTALL_TESTS)
		set(TESTS_EXCLUDE_OPTION "")
	else()
		set(TESTS_EXCLUDE_OPTION EXCLUDE_FROM_ALL)
	endif()

endif()

# configuration scripts (split into multiple files for cleanliness)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(BUILD_HOST_TOOLS) #include only the couple of files needed for building the host tools

	#now enable the languages
	enable_language(C CXX)
	
	include(AmberBuildSystem2ndInit)

	include(CppVectorGen)
	include(ExternalPrograms)

	set(NEEDED_3RDPARTY_TOOLS libm)
	set(REQUIRED_3RDPARTY_TOOLS "")
	set(LINALG_LIBS_REQUIRED FALSE)
	
	include(3rdPartyTools) 
	include(CompilerFlags)
	include(RunnablePrograms)
	
	#CPack setup
	# --------------------------------------------------------------------
	set(PACKAGE_NAME "Amber Host Tools")
	set(PACKAGE_FILENAME "AmberHostTools")
	set(BUNDLE_IDENTIFIER org.ambermd.amber_host_tools)
	set(BUNDLE_SIGNATURE AMHT)
	include(Packaging)
	
	# Tools we need to build: ucpp, rule_parse, nab2c, and utilMakeHelp
	set(AMBER_TOOLS ucpp-1.3 leap)
	
else() #include all configuration segments
	
	#now enable the languages
	enable_language(C CXX Fortran)
	
	include(AmberBuildSystem2ndInit)	
	
	# code for using built-in tools
	# include(nab)
	
	#configuration code
	include(OpenMPConfig)
	include(MPIConfig)
	include(CppVectorGen)
	include(CudaConfig)
	include(PythonInterpreterConfig)
	include(PythonBuildConfig)
	include(ExternalPrograms)
	include(GUIConfig)
	include(PerlConfig)
	include(RunnablePrograms)
	include(AmberUpdater)

    # TCPB
    option(BUILD_TCPB "Enable building TCPB, protobuf library, and TCPB's features in sander." FALSE)
    if(BUILD_TCPB)
	    set(BUILD_PROTOBUF TRUE)
    endif()

	# 3rd party tools
	# --------------------------------------------------------------------

	# These lists control the 3rd party library finder.

	# First, this has all of the tools to search for and potentially use.
	set(NEEDED_3RDPARTY_TOOLS
		blas
		lapack
		arpack 
		ucpp
		c9x-complex
		netcdf
		netcdf-fortran
		protobuf
		fftw
		readline  
		xblas
		lio
		apbs
		pupil
		zlib
		libbz2
		plumed
		libm
		mkl
		boost
		mbx
		kmmd
		tng_io
		nlopt
		libtorch)

	# nccl is only used in pmemd.cuda
	if(CUDA AND NOT AMBERTOOLS_ONLY)
		list(APPEND NEEDED_3RDPARTY_TOOLS nccl)
	endif()

	# some tools are only needed for MPI
	if(MPI)
		if(BUILD_PYTHON AND NOT CROSSCOMPILE)
			list(APPEND NEEDED_3RDPARTY_TOOLS mpi4py)
		endif()
		list(APPEND NEEDED_3RDPARTY_TOOLS pnetcdf)
	endif()

	if(BUILD_PERL)
		list(APPEND NEEDED_3RDPARTY_TOOLS perlmol)
	endif()

	# Next, this has the tools that must be found, internally or externally,
	# for Amber to successfully configure.
	set(REQUIRED_3RDPARTY_TOOLS ucpp boost netcdf netcdf-fortran xblas) 
	if(USE_FFT)
		list(APPEND REQUIRED_3RDPARTY_TOOLS fftw)
	endif()

	if(BUILD_PROTOBUF)
		list(APPEND REQUIRED_3RDPARTY_TOOLS protobuf)
		# Some libraries might be installed at "lib64" instead of "lib", so avoid this from happening
        set(CMAKE_INSTALL_LIBDIR "lib")
	endif()

        if(USE_XTB)
          list(APPEND NEEDED_3RDPARTY_TOOLS xtb)
        endif()
	
        if(USE_DFTBPLUS)
          list(APPEND NEEDED_3RDPARTY_TOOLS dftbplus)
        endif()

	if(USE_DEEPMDKIT)
	  list(APPEND NEEDED_3RDPARTY_TOOLS deepmd-kit)
	endif()
	
	# whether blas and lapack are required
	set(LINALG_LIBS_REQUIRED TRUE)
	
	# This lists the programs that are bundled inside AmberTools/src.
	set(BUNDLED_3RDPARTY_TOOLS blas lapack arpack ucpp c9x-complex kmmd netcdf netcdf-fortran pnetcdf protobuf fftw xblas mpi4py perlmol boost readline tng_io nlopt libtorch)

	# Finally, this selects the libraries to not select as external by default.
	set(SUSPICIOUS_3RDPARTY_TOOLS mkl boost)

	include(3rdPartyTools)
	
	# --------------------------------------------------------------------
	
	include(SanderConfig)
	include(CompilerFlags)
	include(PMEMDCompilerFlags)
	include(ModulePaths)
	include(NetlibConfig)
	
	#tools configuration code
	include(WhichTools)

	# packaging config
	include(cmake-packaging/AmberPackaging.cmake)

	#generate amber.sh
	include(cmake-packaging/PathScript.cmake)

	#install licenses
	include(cmake-packaging/AmberLicenses.cmake)
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
message(STATUS "Configuring subdirs...")

add_subdirectory(AmberTools)

if(NOT BUILD_HOST_TOOLS)
	
	if(NOT AMBERTOOLS_ONLY)
		add_subdirectory(src)
		add_subdirectory(benchmarks)
	endif()
	
	if(NOT AMBER_RELEASE)
		add_subdirectory(examples)
	endif()
	
	add_subdirectory(dat)
	
	# these directories are installed in some configurations
	add_subdirectory(test)
	
	# install compatibility config.h so that tests work
	include(cmake-packaging/LegacyConfigHeader.cmake)
	
	# the doc subdirectory is only present when Amber has been cloned from git.  It does not exist in release builds.
	# Release builds, on the other hand, should have a prebuilt Amber<version>.pdf that we can install
	
	if(AMBER_RELEASE)
		if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/doc/Amber${${PROJECT_NAME}_MAJOR_VERSION}.pdf)
			install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/doc/Amber${${PROJECT_NAME}_MAJOR_VERSION}.pdf DESTINATION ${DOCDIR} COMPONENT Documentation)
		endif()
		
	else()
		include(DocConfig)
		add_subdirectory(doc)
	endif()
endif()

# after configuring all subdirs, install miniconda, since its install rules need to run last
if(BUILD_PYTHON AND DOWNLOAD_MINICONDA)
	install_miniconda()
endif()

# Generate install rule to bundle and fix OS X dependencies
if(NOT $ENV{CONDA_BUILD})
	add_subdirectory(cmake/BundleOSXDependencies)
endif()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

print_build_report()

# print final instructions to user
message(STATUS "Environment resource files are provided to set the proper environment")
message(STATUS "variables to use AMBER and AmberTools. This is required to run any Python")
message(STATUS "programs (like MMPBSA.py, ParmEd, MCPB.py, and pytraj)")
message(STATUS "")

if(TARGET_WINDOWS)
	message(STATUS "You may run the ${CMAKE_INSTALL_PREFIX}/amber.bat file from a CMD shell,")
	message(STATUS "or double-click ${CMAKE_INSTALL_PREFIX}/amber-interactive.bat from Windows")
	message(STATUS "Explorer to start a configured shell for you.")
else()
	message(STATUS "If you use a Bourne shell (e.g., bash, sh, zsh, etc.), source the")
	message(STATUS "${CMAKE_INSTALL_PREFIX}/amber.sh file in your shell. Consider adding the line")
	message(STATUS "  test -f ${CMAKE_INSTALL_PREFIX}/amber.sh && source ${CMAKE_INSTALL_PREFIX}/amber.sh")
	message(STATUS "to your startup file (e.g., ~/.bashrc)")
	message(STATUS "")
	message(STATUS "If you use a C shell (e.g., csh, tcsh), source the")
	message(STATUS "${CMAKE_INSTALL_PREFIX}/amber.csh file in your shell. Consider adding the line")
	message(STATUS "  test -f ${CMAKE_INSTALL_PREFIX}/amber.csh && source ${CMAKE_INSTALL_PREFIX}/amber.csh")
	message(STATUS "to your startup file (e.g., ~/.cshrc)")
	message(STATUS "")
endif()

if(TARGET_OSX OR HOST_OSX)
	message(STATUS "NOTE: MacOS users might need to source it in ~/.bash_profile file")
	message("")
endif()

# be sure to print "mingw32-make" on Windows etc.
get_filename_component(MAKE_PROG_NAME ${CMAKE_MAKE_PROGRAM} NAME)

message(STATUS "Amber will be installed to ${CMAKE_INSTALL_PREFIX}")

#if(CMAKE_GENERATOR MATCHES "Ninja" OR CMAKE_GENERATOR MATCHES "Makefiles")
#	message(STATUS "The next step is to run \"${MAKE_PROG_NAME} install\"")
#endif()
