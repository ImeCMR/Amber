#!/bin/sh
# Touch relevant files to avoid out-of-date issues in Git.
# Run this script in $AMBERHOME/AmberTools/src after switching branches.
# This file is NOT required for the release.
FixTimestamps() {
    sleep 1
    # aclocal-generated aclocal.m4 depends on locally-installed
    # '.m4' macro files, as well as on 'configure.ac'
    touch aclocal.m4
    sleep 1
    # autoconf-generated configure depends on aclocal.m4 and on
    # configure.ac
    touch configure
    # so does autoheader-generated config.h.in
    touch config.h.in
    sleep 1
    touch config.status
    # and all the automake-generated Makefile.in files
    touch `find . -name Makefile.in -print`
    sleep 1
    touch Makefile
}
echo "Warning: NetCDF build timestamps need to be updated. This should"
echo "         ONLY happen during development after switching branches"
echo "         in GIT. Fixing timestamps (this will take about 7 seconds)."
cd netcdf-4.3.0
FixTimestamps
touch libsrc/netcdf.3
cd ../netcdf-fortran-4.2
FixTimestamps
touch man4/stamp-1 man4/version-f90.texi man4/stamp-vti man4/version-f77.texi
touch man4/*.info*
cd ..
