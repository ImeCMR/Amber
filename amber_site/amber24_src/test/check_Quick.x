#!/bin/csh -f
#checks if the file quick exists
#if it does, it is assumed that we have a valid Quick installation

which quick >& /dev/null
if( $status ) then
    # we didn't exit, so Quick seems not installed...
    echo 'Quick not installed - Skipping Test...'
    echo ''
    exit(1)
endif

exit(0)

