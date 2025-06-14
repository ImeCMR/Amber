include(m4/cblas.m4)dnl
dnl M4 file for generating top level Makefile
dnl
dnl To generate a C-source only version (without reference to M4 stuff)
dnl define the symbol C_only (run m4 with -D C_only).
dnl
dnl To generate a version without plain BLAS routines (i.e., only generate
dnl mixed and extended precision version), define the symbol no_plain_blas
dnl (run m4 with -D no_plain_blas).
dnl
dnl To generate version appropriate for LAPACK inclusion, define the 
dnl symbol LAPACK (run m4 with -D LAPACK).
dnl
dnl To add an entry, add the name to the macro LIBNAMES in cblas.m4.
dnl If there is any dependencies between test libraries, add it
dnl in the macro TESTLIB_DEPENDENCY below.
dnl
dnl
define(`TESTLIB_DEPENDENCY', `
# Test library dependencies
sum-test-lib: dot-test-lib
axpby-test-lib: dot-test-lib
waxpby-test-lib: dot-test-lib
gemv-test-lib: dot-test-lib
ge_sum_mv-test-lib: dot-test-lib gemv-test-lib symv-test-lib gemm-test-lib
gbmv-test-lib: dot-test-lib
symv-test-lib: dot-test-lib gemm-test-lib
spmv-test-lib: dot-test-lib symv-test-lib
sbmv-test-lib: dot-test-lib symv-test-lib
hemv-test-lib: dot-test-lib symv-test-lib
hbmv-test-lib: dot-test-lib symv-test-lib hemv-test-lib sbmv-test-lib
hpmv-test-lib: dot-test-lib symv-test-lib hemv-test-lib
trmv-test-lib: dot-test-lib
tpmv-test-lib: dot-test-lib
trsv-test-lib: dot-test-lib trmv-test-lib
tbsv-test-lib: dot-test-lib gbmv-test-lib trsv-test-lib
gemm-test-lib: dot-test-lib gemv-test-lib
symm-test-lib: dot-test-lib gemv-test-lib symv-test-lib
hemm-test-lib: dot-test-lib gemv-test-lib hemv-test-lib symv-test-lib symm-test-lib
gemv2-test-lib: dot2-test-lib gemv-test-lib
dot2-test-lib: dot-test-lib
')dnl
dnl
dnl
ifdef(`LAPACK', `define(`C_only')define(`no_plain_blas')')dnl
define(`if_m4', `ifdef(`C_only', `$2', `$1')')dnl
define(`if_lapack', `ifdef(`LAPACK', `$1', `$2')')dnl
dnl
dnl
changecom()dnl disable comments starting with #
dnl
dnl
define(`MAKEFILE_ENTRY', `
# $1 stuff

$1: $1-test

if_m4(`
$1-source:
	cd $(M4_DIR)/$1 && $(MAKE)
')dnl

$1-lib: common-lib`'if_m4(` $1-source')
	cd $(SRC_DIR)/$1 && $(MAKE)

$1-lib-amb: common-lib`'if_m4(` $1-source')
	cd $(SRC_DIR)/$1 && $(MAKE) all-amb
if_m4(`
$1-test-source:
	cd $(M4_DIR)/test-$1 && $(MAKE)
')dnl

$1-test-lib: if_m4(`$1-test-source ')lib common-test-lib
	cd $(TEST_DIR)/test-$1 && $(MAKE) do_test_$1 

$1-test: $1-test-lib
	cd $(TEST_DIR)/test-$1 && $(MAKE)

$1-clean:
	cd $(SRC_DIR)/$1 && $(MAKE) clean
	cd $(TEST_DIR)/test-$1 && $(MAKE) clean

if_m4(`
$1-source-clean:
	cd $(M4_DIR)/$1 && $(MAKE) source-clean
	cd $(M4_DIR)/test-$1 && $(MAKE) source-clean
')dnl
')dnl
dnl
dnl
#
# This file autogenerated from Makefile.m4 by configure. Do not edit.  
#
# To add a subroutine to be included in libxblas-amb.a:
# 1) Add <subroutine> to the list LIBNAMES_AMB in m4/cblas.m4
# 2) For each new subroutine added, edit m4/<subroutine>/<subroutine>-common.m4
#    to define a new list call <subroutine>_AMB_ARGS.  For an example,
#    see m4/gemv/gemv-common.m4
# 3) Run 'make clean'
# 4) Run 'make makefiles' to generate the new Makefiles in src with
#    the correct Amber targets.
# 5) Run 'make lib-amb' to build libxblas-amb.a
# 6) If you get an error like
#
#    make[1]: *** No rule to make target `BLAS_GEMV2_AMB_ARGSgemv2__.o', needed by `all-amb'.  Stop.
#
#    Redo steps 3), 4), and 5).
#
# Top-level makefile for XBLAS.
#
# To generate all if_m4(`sources, ')libraries`'if_m4(`,') and test results, just type `make'
#
if_m4(`# To generate just the sources, type `make sources'
#
')dnl
# To generate if_m4(`the sources and ')the library, type `make lib'
#
if_m4(`# To generate the test sources, type `make test-sources'
#
')dnl
# To generate all if_m4(`the sources, ')the library, and the test library, type `make test-lib'
#
if_m4(`# To generate all the makefiles but this one, type `make makefiles'
#
')dnl
if_m4(`# To generate this Makefile, type `make Makefile'
#
')dnl
# To clean out all the object files, type `make clean'
#
if_m4(`# To clean out all the object files and source files, type `make source-clean'
#
')dnl
`include' make.conf
`include' $(MAKEINC)

LIB =FOREACH(`LIBNAMES', ` arg-lib')
if_m4(`SOURCES =FOREACH(`LIBNAMES', ` arg-source')
')dnl
LIB_AMB =FOREACH(`LIBNAMES_AMB', ` arg-lib-amb')
if_m4(`SOURCES =FOREACH(`LIBNAMES_AMB', ` arg-source')
')dnl
TESTS =FOREACH(`LIBNAMES', ` arg-test')
if_m4(`TEST_SOURCES =FOREACH(`LIBNAMES', ` arg-test-source') dot2-test-source
')dnl
TEST_LIB =FOREACH(`LIBNAMES', ` arg-test-lib')
CLEAN =FOREACH(`LIBNAMES', ` arg-clean') dot2-clean
if_m4(`SOURCE_CLEAN =FOREACH(`LIBNAMES', ` arg-source-clean') dot2-source-clean
')dnl

SRC_DIR = src
if_m4(`M4_DIR = m4
')dnl
TEST_DIR = testing

all: tests
if_m4(`
sources: $(SOURCES) header

header:
	cd m4 && $(MAKE) header

test-sources: $(TEST_SOURCES)
')dnl

$(TEST_LIB): lib

test-lib: if_m4(`test-sources ')lib $(TEST_LIB)

$(TESTS): test-lib

tests: test-lib $(TESTS)
	rm -f testall.result testall.summary
FOREACH(`LIBNAMES', 
`	cat $(TEST_DIR)/test-arg/arg.results >> testall.result
')dnl
	grep 'FAIL/TOTAL' testall.result >testall.summary
	cat testall.summary

common-lib:
	cd $(SRC_DIR)/common && $(MAKE)

common-test-lib:
	cd $(TEST_DIR)/common && $(MAKE)

if_m4(`
makefiles: Makefile m4/Makefile src/Makefile
	cd m4 && $(MAKE) makefiles
	cd src && $(MAKE) makefiles

src/Makefile: src/Makefile.m4
	cd src && $(M4) $(M4_OPTS) Makefile.m4 >Makefile

m4/Makefile: m4/Makefile.m4 m4/cblas.m4
	cd m4 && $(M4) $(M4_OPTS) Makefile.m4 >Makefile
')dnl
lib: if_m4(`sources ')$(LIB)
	cd $(SRC_DIR)/common && $(MAKE) lib
FOREACH(`LIBNAMES', 
`	cd $(SRC_DIR)/arg && $(MAKE) lib
')dnl

lib-amb: if_m4(`sources ')$(LIB_AMB)
	cd $(SRC_DIR)/common && $(MAKE) lib-amb
FOREACH(`LIBNAMES_AMB', 
`	cd $(SRC_DIR)/arg && $(MAKE) lib-amb
')dnl

if_m4(`
Makefile: Makefile.m4
	$(M4) $(M4_OPTS) Makefile.m4 > Makefile
')dnl
dnl
dnl Custom rules for test-dot2.
dnl
# custom test-dot2 stuff

if_m4(`
dot2-test-source:
	cd $(M4_DIR)/test-dot2 && $(MAKE)
')dnl
dot2-test-lib: lib
	cd $(TEST_DIR)/test-dot2 && $(MAKE) do_test_dot2 

dot2-test: dot2-test-lib
	cd $(TEST_DIR)/test-dot2 && $(MAKE)

dot2-clean:
	cd $(TEST_DIR)/test-dot2 && $(MAKE) clean

if_m4(`
dot2-source-clean:
	cd $(M4_DIR)/test-dot2 && $(MAKE) source-clean
')dnl
dnl
dnl
FOREACH(`LIBNAMES', `MAKEFILE_ENTRY(arg)')dnl

TESTLIB_DEPENDENCY

# Cleaning stuff

clean: $(CLEAN)
	cd $(SRC_DIR)/common && $(MAKE) clean
	cd $(TEST_DIR)/common && $(MAKE) clean
if_m4(`	cd $(M4_DIR) && $(MAKE) clean
')dnl

dist-clean: clean
	rm -f $(XBLASLIB) testall.result testall.summary
if_lapack(`', `	rm -rf autom4te.cache
	rm -f config.log config.status
')dnl

if_m4(`
source-clean: dist-clean
	cd m4 && $(MAKE) source-clean

maintainer-clean: source-clean
	cd m4 && $(MAKE) maintainer-clean
	rm -f m4/Makefile
	cd src && $(MAKE) maintainer-clean
	rm -f src/Makefile
	rm -f Makefile
	rm -f configure

')dnl
.PHONY: $(LIBS) $(TEST_LIB) $(TESTS) $(CLEAN) \
if_m4(`        $(SOURCES) $(TEST_SOURCES) $(SOURCE_CLEAN) \
        sources test-sources source-clean maintainer-clean \
')dnl
        all test-lib tests common-lib common-test-lib lib lib-amb clean dist-clean

