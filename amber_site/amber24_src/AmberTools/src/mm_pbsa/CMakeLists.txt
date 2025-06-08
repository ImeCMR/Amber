project(mm_pbsa)

add_subdirectory(Examples)

set(MAKE_CRD_HG_SOURCES make_crd_hg.F)

#perl scripts to install
set(MM_PBSA_PERL_PROGS mm_pbsa.pl mm_pbsa_statistics.pl)
set(MM_PBSA_PERL_MODULES mm_pbsa_createinput.pm  mm_pbsa_statistics.pm
	mm_pbsa_calceneent.pm mm_pbsa_global.pm mm_pbsa_util.pm
	mm_pbsa_createcoords.pm  mm_pbsa_readinput.pm ForkManager.pm)

#------------------------------------------------------------------------------------------

# I'm assuming that the mm_pbsa programs should be optimized.

# make_crd_hg executable
add_executable(make_crd_hg ${MAKE_CRD_HG_SOURCES})
set_property(TARGET make_crd_hg PROPERTY COMPILE_OPTIONS ${OPT_FFLAGS})
target_link_libraries(make_crd_hg amber_common)

# mm_pbsa nabnmode
include_directories(../nabc ../sff)
add_executable(mm_pbsa_nabnmode mm_pbsa_nabnmode.c)
target_link_libraries(mm_pbsa_nabnmode nabc sff libpbsa cifparse)

install(TARGETS make_crd_hg mm_pbsa_nabnmode DESTINATION ${BINDIR})
install(PROGRAMS ${MM_PBSA_PERL_PROGS} DESTINATION ${BINDIR})
install(FILES ${MM_PBSA_PERL_MODULES} DESTINATION ${PERL_MODULE_INSTALL_DIR}/mm_pbsa)

