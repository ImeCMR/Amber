#!/bin/bash
set -e
set -u

echo ""
echo "Running classical MD in the NVT ensemble to equilibrate"
echo "the system after changing water models. In principle, one"
echo "should also requilibrate the system density, and the "
echo "simulations should be run longer than what this example"
echo "demonstrates."
echo ""

sander -O -p qspcfw.parm7 -i nvt.mdin -o nvt.mdout \
       -c qspcfw.rst7 -r nvt.rst7 -x nvt.nc -inf nvt.mdinfo


echo ""
echo "Brief run of classical MD in the NVE ensemble using sander"
echo ""

sander -O -p qspcfw.parm7 -i nve.mdin -o nve.mdout \
       -c nvt.rst7 -r nve.rst7 -x nve.nc -inf nve.mdinfo



