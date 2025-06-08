#!/bin/bash

echo ""
echo "Running parmutils-rst2ipixyz.py to transform nve.rst7 to nve.xyz"
echo "and create a template i-pi input file called nve.xml."
echo ""

python3 parmutils-rst2ipixyz.py \
	-p qspcfw.parm7 -c nve.rst7 \
	-d system.disang --nve

echo ""
echo "Setting the timestep to 0.5 fs/step because we are using a"
echo "flexible water model"
echo ""

# 8000 steps * 0.5 fs/step = 4 ps

sed -i -e "s|>.*</total_steps>|>8000</total_steps>|" \
    -e "s|>.*</timestep>|>0.5</timestep>|" \
    -e "s|stride='4'|stride='20'|" nve.xml

