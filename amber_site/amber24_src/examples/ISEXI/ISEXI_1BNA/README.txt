RUNNING AN IMPLICIT SOLVENT EXPLICIT ION (ISEXI) SIMULATION OF A SHORT DNA FRAGMENT.

CREDITS:  Yeyue Xiong (code development and testing)  Yegor Kolesnikov (parameter development and performance testing in MD simulations), Alexey Onufriev (research design and project supervision). 

In this simulation we will use the "standard" B-DNA 12 base pair structure 1BNA, neutralized with 22 Na+ ions. 14 NaCl pairs also added to the system to create the ionic atmosphere corresponding to 0.145M [NaCl]. SLTCAP tool (J. Chem. Theory Comput. 14, 4, 1823-1827 (2018)) was used to estimate the numbers of ions. The simulation volume is defined by a spherical potential that prevents the ions from drifting beyond 40 angstroms from the DNA defined as the center of mass of 10 atoms closest to the geometric center of the whole molecule. We are using AMBER's nmropt=1, and its associated parameters, to set the restraints.


1. Make a directory with the files listed below:
	1. {min, nvt, equi, prod}.in
	2. simulation.sh
	3. dna.top dna.crd
	4. disang_NaCl.txt

2. Run the MD simulation on a GPU card. 

2a. Choose the graphics card using command
	export CUDA_VISIBLE_DIVECES=n
Where n is a number of the card used

2b. launch the simulation using BASH script
	./simulation.sh


COMMANDS IN simulation.sh:

2c. pmemd.cuda -O -i min.in -o min.out -p dna.top -c dna.crd -r min.ncrst -inf min.mdinfo -ref dna.crd
	–---
	initially minimize the system
	–---

2d. pmemd.cuda -O -i heat.in -o heat.out -p dna.top -c min.ncrst -r heat.ncrst -x heat.trj -inf heat.mdinfo -ref dna.crd
	–---
	heats the system to 300K
	–---

2e. pmemd.cuda -O -i equil.in -o equil.out -p dna.top -c heat.ncrst -r equil.ncrst -x equil.trj -inf equil.mdinfo
	–---
	equilibrates the system
	–---

2f. pmemd.cuda -O -i prod.in -o prod.out -p dna.top -c equil.ncrst -r prod.ncrst -x prod.trj -inf prod.mdinfo
	–---
	launch 10 nanosecond long production simulation
	–---

2g. ambpdb -p dna.top -c prod.ncrst > prod.pdb
	–---
	generate pdb file with result of the simulation
	–---

3. RESULTS. 

The results of the 10 nanosecond production run are in prod.{out,trj,ncrst} files. The final system snapshot is prod.pdb.

prod.ncrst file is a restart file for further simulation if necessary
prod.trj file contain coordinates of atoms in every snapshot
prod.out file contain energies of the system in every snapshot

