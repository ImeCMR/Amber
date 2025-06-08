#!/bin/bash
set -e
set -u

echo ""
echo "Starting i-pi..."
i-pi nve.xml &> nve.stdout &


echo ""
echo "Waiting for initialization..."
sleep 10

echo ""
echo "Starting sander..."
mpirun -n 4 sander.MPI -host 127.0.0.1 -port 31415 \
       -O -p qspcfw.parm7 -c nve.rst7 \
       -i ipi.mdin -o ipi.mdout -r ipi.rst7 \
       -x ipi.nc -inf ipi.mdinfo &

wait

echo ""
echo "Finished!"
echo ""
echo "To see the total energy versus time, try running:"
echo "   xmgrace -block nve.out -bxy 2:4"
echo "The x-axis is time (ps) and the y-axis is energy (cal/mol)"
echo ""


