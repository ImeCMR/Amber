set(NMODE_SOURCES nmode.F90 alloc.F90 angl.F90 assert.F90 
	bond.F90 bresna.F90 cdiv.F90 chck.F90 chdir.F90 chlamba.F90 
	debug.F90 decomp.F90 difang.F90 difbon.F90 dihdup.F90 dihed.F90 
	dirch.F90 dsort.F90 
	eigout.F90 ephi.F90 eprnt.F90 errstp.F90 
	findl.F90 forces.F90 
	getcor.F90 getf.F90 getm.F90 getref.F90 grdmax.F90 
	level.F90 lmdout.F90 lmode.F90 
	makel.F90 movecm.F90 mweight.F90 mweit.F90 
	new_stuff.F90 
	nmdfil.F90 nonbon.F90 nonbon1.F90 norm.F90 norml.F90 nrch.F90 
	oseen.F90 pairr.F90 pairw.F90 percon.F90 pload.F90 putvar.F90 
	rdinp.F90 rdparm.F90 remrot.F90 resdis.F90 
	savec.F90 savei.F90  savet.F90 setang.F90 setatm.F90 setbon.F90 
	setcor.F90 setdih.F90 setgam.F90 setgr.F90 setvar.F90 spsort.F90 
	thermo.F90 timit.F90 traco.F90  tsearch.F90 tstate.F90 
	xconst.F90 zerout.F90 zxcgr.F90 mexit.F90)
	
set(NMANAL_SOURCES mdread.F90 nmanal.F90 rdvect.F90 proj.F90 stigrp.F90 setgrp.F90
    group.F90 grpprj.F90 rdparm2.F90 rmsf.F90 excise.F90 dotref.F90
    prtma.F90 reori.F90 movecm.F90 bmat.F90 belly.F90 quasip.F90
    fluct.F90 nmdfil.F90 corf.F90 mexit.F90 lmdin.F90 
    corfl.F90 corffft.F90 corffft2.F90 ../sander/pubfft.F90 fileio.F90 rgroup2.F90)
        
set(LMANAL_SOURCES lmanal.F90 distr.F90 lmfil.F90 mexit.F90 count.F90 
    fluct.F90 lmdin.F90 dsort.F90 getref.F90 fileio.F90 
    corfl.F90 corfv.F90 rdparml.F90)
    
set(QUASIH_SOURCES quasih.F90 fileio.F90 thermo.F90 mexit.F90 decomp.F90 assert.F90)

set(CALINK_SOURCES calink.F90 fileio.F90)

set(VECS_TO_BINPOS_SOURCES vecs_to_binpos.F90 fileio.F90)

#source files which should be compiled with optimizations.  All other soures are unoptimized.
set(OPTIMIZED_SOURCES calink.F90 rdparm.F90 rdparm2.F90 nonbon.F90 nonbon1.F90 bond.F90 angl.F90 bresna.F90 resdis.F90 ../sander/pubfft.F90)

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

add_compile_options(${NO_OPT_FFLAGS})

set_property(SOURCE ${OPTIMIZED_SOURCES} PROPERTY COMPILE_OPTIONS ${OPT_FFLAGS})

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------

add_executable(nmode ${NMODE_SOURCES})
config_module_dirs(nmode ${NMODE_MOD_DIR})

# add_executable(nmanal ${NMANAL_SOURCES})
# config_module_dirs(nmanal ${NMODE_MOD_DIR})

# add_executable(lmanal ${LMANAL_SOURCES})
# config_module_dirs(lmanal ${NMODE_MOD_DIR})

# add_executable(quasih ${QUASIH_SOURCES})

# add_executable(calink ${CALINK_SOURCES})

# add_executable(vecs_to_binpos ${VECS_TO_BINPOS_SOURCES})

targets_link_libraries(nmode LIBRARIES netlib amber_common)

install(TARGETS nmode DESTINATION ${BINDIR})
