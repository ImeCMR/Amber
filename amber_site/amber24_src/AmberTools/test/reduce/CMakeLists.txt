# CMake test file for reduce

add_dacdif_test(reduce-test
		COMMAND $<TARGET_FILE:reduce> -build 1ubq.pdb
		OUTPUT_FILES @stdout@ @stderr@
		SAVE_FILES 1ubqH.pdb.save messages.save)