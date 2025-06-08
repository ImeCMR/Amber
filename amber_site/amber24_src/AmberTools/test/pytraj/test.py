from __future__ import print_function
import pytraj as pt
from pytraj.testing import get_fn
from pytraj.testing import aa_eq

print(pt)

traj = pt.iterload(*get_fn('rna'))
state = pt.load_cpptraj_state('''
        rms @P first
        radgyr nomax''', traj)

state.run()

try:
    assert traj.n_atoms == 518, 'n_atoms = 518'
    assert traj.n_frames == 3, 'n_frames = 3'
    aa_eq(pt.rmsd(traj, mask='@P'), state.data[1].values)
    aa_eq(pt.radgyr(traj), state.data[2].values)
except AssertionError:
    print('Possible FAILURE')
else:
    print('PASSED')
