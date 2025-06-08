#!/bin/bash
set -e
set -u

echo ""
echo "Running tleap on tleap.in to build qspcfw.parm7 and qspcfw.rst7"
echo "from MOL.pdb, MOL.lib, and MOL.frcmod"
echo ""

tleap -s -f tleap.in

echo ""
echo "Running ChBox to set the periodic box, which was previously"
echo "equilibrated using TIP4P/Ew waters"
echo ""

ChBox -X 34.9384593 -Y 37.1182234 -Z 35.4310664 \
      -al 90. -bt 90. -gm 90. \
      -c qspcfw.rst7 -o tmp.rst7

mv tmp.rst7 qspcfw.rst7



