#!/bin/sh

echo ""
echo "--------------------------------------------------------------------------------"
echo "A required NetCDF file '$1' is missing."
echo "Either re-run 'configure' in \$AMBERHOME or check that your NetCDF distribution"
echo "is correctly installed (it must have support for C and Fortran)."
echo "--------------------------------------------------------------------------------"
echo ""
exit 1
