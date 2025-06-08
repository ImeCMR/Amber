# CMake test file for ambpdb

add_dacdif_test(ambpdb-pdb
		COMMAND $<TARGET_FILE:ambpdb> -p prmtop
		STDIN_FILE restrt
		OUTPUT_FILES @stdout@
		SAVE_FILES restrt.pdb.save)

add_dacdif_test(ambpdb-mol2
		COMMAND $<TARGET_FILE:ambpdb> -p prmtop -mol2
		STDIN_FILE restrt
		OUTPUT_FILES @stdout@
		SAVE_FILES restrt.mol2.save)

add_dacdif_test(ambpdb-pqr
		COMMAND $<TARGET_FILE:ambpdb> -p prmtop -pqr
		STDIN_FILE restrt
		OUTPUT_FILES @stdout@
		SAVE_FILES restrt.pqr.save)