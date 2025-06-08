project(nabc)

set(NABC_SOURCES axis2frame.c	
	bonded_atoms.c 
	database.c  
	deform_ener.c 
    fd_helix.c
	frameutil.c	
	getdist.c	
	jacobi.c	
	match.c		
	matop.c		
	molio.c		
	molsurf.c   
	molutil.c	
	pointutil.c	
	prm.c 
	regexp.c 
	reslib.c	
	select_atoms.c	
    stringutil.c 
	traceback.c
    wc_basepair.c
    wc_complement.c
    wc_helix.c)
	
# --------------------------------------------------------------------

add_compile_options(${OPT_CFLAGS})
include_directories(../sff)

# --------------------------------------------------------------------

add_library(nabc STATIC ${NABC_SOURCES})
install_libraries(nabc)

install(FILES nabc.h DESTINATION ${INCDIR})
