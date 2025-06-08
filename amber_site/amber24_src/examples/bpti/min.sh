# Running minimization for BPTI

cat << eof > min.in
# 200 steps of minimization, generalized Born solvent model
&cntrl
maxcyc=200, imin=1, cut=12.0, igb=1, ntb=0, ntpr=10,
/
eof

sander -i min.in -o 6pti.min1.out -c prmcrd -r 6pti.min1.xyz
/bin/rm min.in

