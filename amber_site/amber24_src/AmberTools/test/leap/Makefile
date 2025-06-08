# Makefile for tleap tests

# The default test target
test: test.all

# ----- Individual test targets ------------------
# Test various types of parameter assignment.
# These are part of the test.params target.
# The OPT variable is used to pass arguments to TestCommon.sh, which
# the run scripts for all tests in this section should be sourcing.
PARAM_TESTS=test.ol3 \
            test.yil \
            test.ff14sb.ol3.bsc0 \
            test.bsc0 \
            test.bsc1 \
            test.ol15 \
            test.pbradii \
            test.rna.ljbb

test.ol3:
	@-cd OL3 && ./Run.OL3.sh $(OPT)

test.yil:
	@-cd Yildirim && ./Run.yildirim.sh $(OPT)
	@-cd Yildirim && ./Run.GACU.sh $(OPT)

test.ff14sb.ol3.bsc0:
	@-cd 4PY5 && ./Run.4PY5.sh $(OPT)

test.bsc0:
	@-cd BSC0 && ./Run.BSC0.sh $(OPT)

test.bsc1:
	@-cd BSC1 && ./Run.BSC1.sh $(OPT)

test.ol15:
	@-cd OL15 && ./Run.OL15.sh $(OPT)

test.pbradii:
	@-cd pbradii && ./Run.pbradii.sh $(OPT)

test.rna.ljbb:
	@-cd RNA.LJbb && ./Run.RNA.LJbb.sh $(OPT)

# ----- Targets for collections of tests ---------
fast:
	@-./Run.tleap
	@-./Run.cmds
	@-./Run.mbondi3
	@-./Run.loadpdb

test.leap: fast
	@-./Run.parm7
	@-./Run.addions
	@-cd glycam/06j && ./Run.glycam
	@-cd glycam/06EPb && ./Run.glycam
	@-cd glycam/06j_10 && ./Run.glycam
	@-cd glycam/06j_12SB && ./Run.glycam
	@-cd BondAugmentations/ && ./Run.1p7e

test.params: $(PARAM_TESTS)

test.all:
	$(MAKE) test.leap test.params

# ----- Other targets ----------------------------
clean:
	@echo "Cleaning all test directories:"
	(find . -name "*.dif" -o -name 'leap.log' -o -name 'leap.out' | \
	 while read dif ;\
	 do \
	     rm -f $$dif ;\
	 done ;\
	)
	-/bin/rm -f x.lib lib94.out ions94.log ff94_prm.out ff94_prm.log ff94.log
	$(MAKE) test.params OPT=clean

check:
	(find . -name "*.dif" ! -size 0c -print | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	)

