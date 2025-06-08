SHELL=/bin/sh
AMBER=Amber24
DIFF_FILE_EXTENSION=dif
OUTPUT_FILE_EXTENSION=out

bench.sander: 
	@echo "Starting benchmarking of ${AMBER} at `date`".
	cd factor_ix; ./bench.factor_ix
#	cd hb; ./bench.hb
	cd jac; ./bench.jac
	cd gb_mb; ./bench.gb_mb
	cd gb_alp; ./bench.gb_alp
	cd trx; ./bench.trx
	@echo ""
	@echo "Benchmarking of ${AMBER} is complete at `date`."
	@echo ""

#export TESTsander='pmemd'  ; cd gb_alp; ./bench.gb_alp
#PMEMD 10 does not support gbsa > 0
#export TESTsander='pmemd'  ; cd trx; ./bench.trx
#  unknown flag: -idip
bench.pmemd: 
	@echo "Starting benchmarking of ${AMBER} at `date`".
	@echo "TESTsander is 'pmemd'"
	export TESTsander='pmemd'  ; cd factor_ix; ./bench.factor_ix
#	export TESTsander='pmemd'  ; cd hb; ./bench.hb
	export TESTsander='pmemd'  ; cd jac; ./bench.jac
	export TESTsander='pmemd'  ; cd gb_mb; ./bench.gb_mb
	@echo ""
	@echo "Benchmarking of ${AMBER} is complete at `date`."
	@echo ""

bench: 
	@echo "Starting benchmarking of ${AMBER} at `date`".
	@echo "TESTsander is '${TESTsander}'"
	@(find . -perm -100 -name 'bench.*' ! -name '*pmemd*' ! -name '*prowat*' -print | \
    sort | \
	while read bench ;\
	do \
	    echo $$bench ;\
	    dir=`dirname $$bench` ;\
	    file=`basename $$bench` ;\
	    cd $$dir ;\
	    ./$$file ;\
	    cd .. ;\
	    echo ;\
	done ;\
	)
	@echo ""
	@echo "Benchmarking of ${AMBER} is complete at `date`."
	@echo ""

clean:
	find . -name "*.$(DIFF_FILE_EXTENSION)" -exec rm {} \;
	find . -name "*.$(OUTPUT_FILE_EXTENSION)" -exec rm {} \;

check:
	@(find . -name "*.$(DIFF_FILE_EXTENSION)" ! -size 0c -print | sort | \
	while read dif ;\
	do \
	    echo $$dif ;\
	    cat $$dif ;\
	    echo ;\
	done ;\
	success_token='TIMINGS' ;\
	tail_token='Total time' ;\
	find . -name "*.$(OUTPUT_FILE_EXTENSION)" ! -size 0c -print | sort | \
	while read out ;\
	do \
	    echo $$out ;\
	    if grep "$$success_token" $$out > /dev/null ; then \
	        cat $$out | grep "$$tail_token" ;\
	    else \
	        echo "FAILURE in $$out !" ;\
	    fi ;\
	done ;\
	)

summary:
	@(success_token='TIMINGS' ;\
	tail_token='A V E R A G E S   O V E R' ;\
	find . -name "*.$(OUTPUT_FILE_EXTENSION)" ! -size 0c -print | sort | \
	while read out ;\
	do \
	    echo $$out ;\
	    if grep "$$success_token" $$out > /dev/null ; then \
	        cat $$out | sed -n "/$$tail_token/,/Nonsetup wallclock/p" ;\
	    else \
	        echo ;\
	        echo "FAILURE in $$out !  Tail of output file:" ;\
	        tail $$out ;\
	    fi ;\
	    echo ;\
	done ;\
	)

