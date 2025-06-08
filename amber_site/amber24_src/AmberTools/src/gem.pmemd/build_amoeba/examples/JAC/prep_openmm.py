#!/usr/bin/env python
from __future__ import division, print_function

import simtk.openmm as mm
import simtk.unit as u
import simtk.openmm.app as app

# Parse the PDB file
print('Parsing the PDB file...')
pdb = app.PDBFile('jac.pdb')

# Declare the FF we want to use... in this case amoeba13.xml
print('Loading the force field...')
ff = app.ForceField('amoeba2013.xml')

# Create the System object with it
print('Creating the System...')
system = ff.createSystem(pdb.topology, nonbondedMethod=app.PME,
                         nonbondedCutoff=7*u.angstroms, rigidWater=False,
                         constraints=None)

# Now scan through our forces and change the cutoff for the van der Waals force
# to 9 angstroms, and change our dipole convergence to 1e-6
print('Adjusting the vdW cutoff to 9 Angstroms...')
for force in system.getForces():
    if isinstance(force, mm.AmoebaVdwForce):
        force.setCutoff(9*u.angstroms)
    elif isinstance(force, mm.AmoebaMultipoleForce):
        force.setMutualInducedTargetEpsilon(1e-6)

# Now we are done creating our system. Let's serialize it by writing an XML file
print('Serializing the System...')
with open('jac.system.xml', 'w') as f:
    f.write(mm.XmlSerializer.serialize(system))

# Finish progress report
print('Done.')
