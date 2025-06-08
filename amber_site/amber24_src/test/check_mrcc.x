#!/bin/csh -f
#checks if the file dmrcc exists and is executable
#if it does, it is assumed that we have a valid MRCC installation

which mrcc >& /dev/null
if( $status ) then
    # we didn't exit, so MRCC seems not installed...
    echo 'MRCC not installed - Skipping Test...'
    echo ''
    exit(1)
endif

exit(0)
