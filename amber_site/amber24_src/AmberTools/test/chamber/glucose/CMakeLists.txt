# CMake test file for chamber
add_dacdif_test(chamber-glucose
		COMMAND $<TARGET_FILE:chamber>
		-top top_cbh1_cellulose.rtf
		-param par_cbh1_cellulose.prm
		-psf gluc.psf -nocmap
		-crd glucmin.pdb
		-p ${CMAKE_CURRENT_BINARY_DIR}/prmtop
		-inpcrd ${CMAKE_CURRENT_BINARY_DIR}/inpcrd
		DACDIF_ARGS -r 10E-15
		OUTPUT_FILES @stdout@ prmtop inpcrd
		SAVE_FILES chamber.out.save prmtop.save inpcrd.save)