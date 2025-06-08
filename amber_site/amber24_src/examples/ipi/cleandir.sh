#!/bin/bash
set -e
set -u

for f in ipi.mdout leap.log nve.mdinfo nve.mdout \
		   nve.nc nve.restart nve.rst7 nve.stdout \
		   nve.xcrd nve.xc.xyz nve.xml nve.xyz \
		   nvt.mdinfo nvt.mdout nvt.nc nvt.rst7 \
		   qspcfw.parm7 qspcfw.rst7 RESTART \
		   system.disang.plumed.dumpave \
		   system.dumpave system.plumed; do
    
    if [ -e "${f}" ]; then
	rm "${f}"
    fi
    
done

