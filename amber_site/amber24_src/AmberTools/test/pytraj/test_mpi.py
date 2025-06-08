import pytraj as pt
from pytraj.testing import get_fn
from pytraj.testing import aa_eq
from mpi4py import MPI

comm = MPI.COMM_WORLD

traj = pt.iterload(*get_fn('tz2'))

parallel_data = pt.pmap_mpi(pt.rmsd, traj, ref=-3, mask='@CA')

if comm.rank == 0:
    serial_data = pt.rmsd(traj, ref=-3, mask='@CA') 
    try:
        aa_eq(parallel_data['RMSD_00001'], serial_data)
    except AssertionError:
        print('Possible FAILURE')
    else:
        print('PASSED')
