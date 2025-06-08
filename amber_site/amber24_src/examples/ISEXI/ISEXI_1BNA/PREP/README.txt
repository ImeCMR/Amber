SETTING UP AN IMPLICIT SOLVENT EXPLICIT ION (ISEXI) SIMULATION OF A SHORT DNA FRAGMENT.

1. Make a directory with the files listed below:
	1. tleap.script
	2. 1bna.pdb
	3. disang.py
	4. prep.sh

2. PREPARE THE INPUT FILES. 

2a. Use prep.sh script, which contains the following commands:

2b tleap -f tleap.script
	–---
	generates topology file and initial coordinates file
	–---

2c python disang_1BNA_NaCl.py
	–---
	generates distance restraints for the ions
	–---

3. RESULTS. 

The results of the script are the files dna.{top,crd} and disang_NaCl.txt.

dna.top – topology file of the system, describes the atoms in the system and their interactions
dna.crd file contain the initial coordinates of atoms in the system
disang_NaCl.txt file contain definitions of the restraining potential forces defining simulation volume

