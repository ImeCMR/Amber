#!/bin/bash

# Test Bergonzo & Cheatham RNA phosphate b.b. params

. ../TestCommon.sh

CleanFiles leap.in leap.out leap.log test.04.GACC.parm7 test.04.GACC.rst7

LEAPRC=leaprc.RNA.LJbb

printf "\nTest Bergonzo & Cheatham phosphate b.b. parameters for RNA.\n\n"
cat > leap.in <<EOF
source $LEAPRC
mol = sequence{G5 A C C3}
set default PBradii mbondi3

saveamberparm mol test.04.GACC.parm7 test.04.GACC.rst7
quit
EOF
 
RunTleap test.04.GACC.parm7

## Test that loading a protein leaprc does not affect Yildirim params
#cat > leap.in <<EOF
#source $LEAPRC
#source leaprc.protein.ff14SB
#m = loadpdb rGACC.pdb
#saveamberparm m rGACC.Yil.parm7 rGACC.Yil.rst7
#quit
#EOF
#RunTleap rGACC.Yil.parm7
$DACDIF test.04.GACC.parm7.save test.04.GACC.parm7
$DACDIF test.04.GACC.rst7.save test.04.GACC.rst7

CleanFiles leap.in leap.out leap.log

exit 0
