#!/bin/bash

# Check the existence of a test version of the mdgx program
MDGX="${AMBERHOME}/bin/mdgx"
if [ -n "${TESTmdgx}" ] ; then
  MDGX=${TESTmdgx}
  echo "MDGX set to ${MDGX}"
fi

# Write a simple input file for plain torsion optimization
cat > mdgx.in << EOF
&files
  -p       ../Torsion/THR.top
  -c       ../Torsion/THR1.crd
  -o       confsamp.out
&end

&configs

  % Operations
  GridPerturb  :2@N :2@CB              {   -0.1     0.2 } Krst 128.0
  GridPerturb  :2@C :2@HG1             {   -1.5     1.5 } Krst 256.0
  GridSample   :2@N :2@CA :2@HA        {  100.0   120.0 } Krst 256.0
  RandomSample :1@C :2@N  :2@CA :2@C   { -160.0  -140.0 } Krst 64.0
  RandomSample :2@N :2@CA :2@CB :2@OG1 {  100.0   140.0 } Krst 32.0 fbhw 10.0

  % Control parameters governing the configuration generation
  count = 8,
  ig = 537618,
  ncyc = 100,
  maxcyc = 1200,

  % Output conditions
  write    'pdb', 'inpcrd',
  outbase  'conf',
  outsuff  'pdb', 'crd',
&end
EOF

# Run the mdgx program
${MDGX} -O -i mdgx.in
grep " E(" confsamp.out > extract.txt
grep "Strain is" confsamp.out >> extract.txt
grep "Restraint binding" confsamp.out >> extract.txt

# Analyze the relevant output
${AMBERHOME}/AmberTools/test/dacdif -r 1.0e-4 extract.txt.save extract.txt
${AMBERHOME}/AmberTools/test/dacdif -r 1.0e-3 conf7.pdb.save conf7.pdb
${AMBERHOME}/AmberTools/test/dacdif -r 1.0e-4 conf5.crd.save conf5.crd

# Try reshuffling
head -21 mdgx.in > tmp
mv tmp mdgx.in
cat >> mdgx.in << EOF

  % Output conditions
  write    'inpcrd',
  outbase  'confr',
  outsuff  'crd',

  % Reshuffle
  nshuffle = 2,
  shuffle 'jackknife',
&end
EOF

# Run the mdgx program
${MDGX} -O -i mdgx.in
grep " E(" confsamp.out > extractr.txt
grep "Strain is" confsamp.out >> extractr.txt
grep "Restraint binding" confsamp.out >> extractr.txt

# Analyze the relevant output
${AMBERHOME}/AmberTools/test/dacdif -r 1.0e-4 extractr.txt.save extractr.txt
${AMBERHOME}/AmberTools/test/dacdif -r 1.0e-4 confr2.crd.save confr2.crd

# Tidy up
/bin/rm -f mdgx.in conf[0-9]*.pdb conf[0-9]*.crd confr[0-9]*.crd 
/bin/rm -f confsamp.out extract.txt
