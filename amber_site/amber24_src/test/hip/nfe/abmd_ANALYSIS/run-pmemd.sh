#!/bin/sh
#TEST-PROGRAM sander
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

. ../common.sh

#
# don't run if there is no netCDF
#

set_NC_UTILS

if test $? -ne 0; then
   exit 0;
fi

set_SANDER pmemd.hip_$1
set_SANDER_CMD 1

JUNK="mdout mdinfo restrt monitor.txt inpcrd prmtop sander.out logfile"

#
# remove the junk
#

/bin/rm -rf ${JUNK} junk.*

#
# prepare files
#

cp -p ../prmtop .
cp -p ../inpcrd.4 ./inpcrd

#
# run SANDER
#

${SANDER_CMD} -O -AllowSmallBox -i mdin > sander.out 2>&1
CheckError $?

if [ "`basename $SANDER`" = "pmemd.hip_DPFP" -o "`basename $SANDER`" = "pmemd.hip_DPFP.MPI" ]; then
  ../../../dacdif -r 1.0e-4 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-5 save.pmemd/monitor.txt monitor.txt
else
  ../../../dacdif -r 1.0e-2 save.pmemd/mdout mdout
  ../../../dacdif -r 1.0e-3 save.pmemd/monitor.txt monitor.txt
fi

#
# preserve the junk on failure
#

save_junk_on_failure ${JUNK}

#
# remove the junk
#

/bin/rm -f ${JUNK} profile_mpi
