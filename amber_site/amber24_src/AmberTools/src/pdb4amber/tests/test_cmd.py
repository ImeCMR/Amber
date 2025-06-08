import os
import sys
import subprocess
import unittest
import parmed as pmd
import numpy as np
from parmed.residue import WATER_NAMES
import pytest
from mock import patch
try:
    from io import StringIO
except ImportError:
    from cStringIO import StringIO

from pdb4amber import pdb4amber
# local
from utils import tempfolder, get_fn, _has_program

pdb_fn = get_fn('4lzt/4lzt_h.pdb')

try:
    # check internet
    pmd.download_PDB('1l2y')
    internet_ok = True
except IOError:
    internet_ok = False


@unittest.skipUnless(internet_ok, 'must have internet connection to rcsb')
def test_write_model():
    fname = get_fn('1l2y.pdb')
    orig_parm = pmd.load_file(fname)
    pdb_out = 'out.pdb'

    # default
    with tempfolder():
        pdb4amber.main([fname, '-o', pdb_out])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (1, 304, 3)

    # model 1
    with tempfolder():
        pdb4amber.main([fname, '-o', pdb_out, '--model', '1'])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (1, 304, 3)
        np.testing.assert_almost_equal(parm.coordinates,
                                       orig_parm.get_coordinates()[0])

    # model 2
    model = 2
    with tempfolder():
        pdb4amber.main([fname, '-o', pdb_out, '--model', str(model)])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (1, 304, 3)
        np.testing.assert_almost_equal(parm.coordinates,
                                       orig_parm.get_coordinates()[model - 1])

    # keep all models
    with tempfolder():
        pdb_out = 'out.pdb'
        pdb4amber.main([fname, '-o', pdb_out, '--model', '-1'])
        parm = pmd.load_file(pdb_out)
        assert parm.get_coordinates().shape == (38, 304, 3)


def test_dry():
    option = '--dry'
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, option]

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        resnames = set(res.name for res in orig_parm.residues)
        assert resnames.intersection(WATER_NAMES)

        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        resnames = set(res.name for res in parm.residues)
        assert not resnames.intersection(WATER_NAMES)

        # water
        water_parm = pmd.load_file('out_water.pdb')
        assert set(res.name for res in water_parm.residues) == {'HOH'}


def test_noh(tmpdir):
    inp_pdb = get_fn('nmr_struc_1.pdb')
    with tmpdir.as_cwd():
        parm_with_h = pmd.load_file(inp_pdb)
        assert len(parm_with_h.atoms) == 304
        assert len([atom for atom in parm_with_h.atoms if atom.atomic_number == 1]) == 150
        pdb4amber.main([inp_pdb, '-y', '-o', 'out.pdb'])
        parm = pmd.load_file('out.pdb')
        assert not [atom for atom in parm.atoms if atom.atomic_number == 1]
        assert len(parm.atoms) == 154


def test_not_write_sslink_conect_record():
    pdb_out = 'out.pdb'
    pdb_fn = get_fn('4lzt/4lzt_h.pdb')
    sslink_name = 'out_sslink'

    # has conect
    with tempfolder():
        command = ['-i', pdb_fn, '-o', pdb_out]
        pdb4amber.main(command)
        with open(pdb_out) as fh:
            assert 'CONECT' in fh.read()

    # no conect
    with tempfolder():
        command = ['-i', pdb_fn, '-o', pdb_out, '--no-conect']
        pdb4amber.main(command)
        with open(pdb_out) as fh:
            assert 'CONECT' not in fh.read()


def test_write_sslink():
    pdb_out = 'out.pdb'
    pdb_fn = get_fn('4lzt/4lzt_h.pdb')
    command = ['-i', pdb_fn, '-o', pdb_out]
    sslink_name = 'out_sslink'
    sslink_pair = [(6, 127), (30, 115), (64, 80), (76, 94)]

    with tempfolder():
        pdb4amber.main(command)
        with open(sslink_name) as fh:
            for index, line in enumerate(fh):
                id0, idx1 = [int(i) for i in line.split()]
                assert (id0, idx1) == sslink_pair[index]


def test_write_conect():
    pdb_out = 'out.pdb'
    pdb_fn = get_fn('1dwc.pdb')
    command = ['-i', pdb_fn, '-o', pdb_out]
    with tempfolder():
        pdb4amber.main(command)
        with open(pdb_out) as fh:
            assert 'ATOM   1196  SG  CYX H 148      55.729  28.382  19.687  1.00 10.07' in fh.read(
            )


def test_constantph():
    option = '--constantph'
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, option]

    with tempfolder():
        # just run to increase code coverage
        # we already test in another file
        pdb4amber.main(command)


def test_no_hydrogen():
    option = '--nohyd'
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, option]

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        atom_names = set(
            atom.name for atom in orig_parm.atoms if atom.atomic_number == 1)
        assert atom_names

        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(
            atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert not atom_names


def test_prot_only(tmpdir):
    option = '--pro'
    pdb_in = get_fn('3orn.pdb')
    pdb_out = 'out.pdb'
    command = ['-i', pdb_in, '-o', pdb_out, option]

    with tmpdir.as_cwd():
        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        res_names = set(res.name for res in parm.residues)
        assert res_names == {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP',
            'TYR', 'VAL'
        }


def test_amber_compatible_residues(tmpdir):
    option = '--amber-compatible-residues'
    pdb_in = get_fn('3orn.pdb')
    pdb_out = 'out.pdb'
    command = ['-i', pdb_in, '-o', pdb_out, option]

    with tmpdir.as_cwd():
        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        res_names = set(res.name for res in parm.residues)
        assert res_names == {
            'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS',
            'HOH', 'ILE', 'LEU', 'LYS', 'MET', 'MG', 'PHE', 'PRO', 'SER',
            'THR', 'TRP', 'TYR', 'VAL'
        }


def test_writing_renum():
    pdb_fn = get_fn('2igd/2igd_2_residues.pdb')
    pdb_out = 'out.pdb'
    command = [pdb_fn]
    expected_lines = """
MET     1    MET     1
PRO     3    PRO     2
    """.strip().split('\n')

    with tempfolder():
        pdb4amber.main(command)
        with open('stdout_renum.txt') as fh:
            output_lines = [line.strip() for line in fh.readlines()]
            assert expected_lines == output_lines


def test_reduce_with_pdb_input():
    option = '--reduce'
    pdb_fn = get_fn('2igd/2igd.pdb')
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, option]

    with tempfolder():
        orig_parm = pmd.load_file(pdb_fn)
        atom_names = set(
            atom.name for atom in orig_parm.atoms if atom.atomic_number == 1)
        assert not atom_names

        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(
            atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert atom_names


def test_strip_atoms():
    pdb_fn = get_fn('2igd/2igd.cif')
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, '--strip', ':3-500']
    with tempfolder():
        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        assert len(parm.residues) == 2


def test_reduce_with_cif_input():
    option = '--reduce'
    pdb_fn = get_fn('2igd/2igd.cif')
    pdb_out = 'out.pdb'
    command = ['-i', pdb_fn, '-o', pdb_out, option]

    with tempfolder():
        pdb4amber.main(command)
        parm = pmd.load_file(pdb_out)
        atom_names = set(
            atom.name for atom in parm.atoms if atom.atomic_number == 1)
        assert atom_names


def test_stdin_stdout():
    ''' e.g: cat my.pdb | pdb4amber '''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['cat', pdb_fn, '|', 'pdb4amber']

    with tempfolder():
        # use shell=True since check_output return exit 1 with |
        # not sure why.
        output = subprocess.check_output(
            ' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 574


@unittest.skipUnless(internet_ok, 'internet')
def test_fetch_pdbid(tmpdir):
    ''' e.g: pdb4amber 1l2y --pdbid '''
    pdb_fn = '1l2y'

    with tmpdir.as_cwd():
        with patch('parmed.download_PDB') as mock_download:

            def effect(*args):
                return pmd.load_file(get_fn('1l2y.pdb'))

            mock_download.side_effect = effect
            pdb4amber.main(['1l2y', '--pdbid', '-o', 'out.pdb'])
            assert mock_download.called


def test_simplest_command_pdb4amber_mypdb():
    # pdb4amber my.pdb
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['pdb4amber', pdb_fn]

    with tempfolder():
        output = subprocess.check_output(
            ' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 574


@pytest.mark.xfail
def test_simplest_command_pdb4amber():
    command = ['pdb4amber']
    with tempfolder():
        output = subprocess.check_output(command).decode()
        assert 'usage: pdb4amber' in output


def test_stdin_stdout_with_reduce():
    ''' e.g: cat my.pdb | pdb4amber --reduce '''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['cat', pdb_fn, '|', 'pdb4amber', '--reduce']

    with tempfolder():
        # use shell=True since check_output return exit 1 with |
        # not sure why.
        output = subprocess.check_output(
            ' '.join(command), shell=True).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        assert len(parm.atoms) == 1033


def test_write_other_formats_like_mol2():
    # mol2
    pdb_out = 'out.mol2'
    command = ['pdb4amber', '-i', pdb_fn, '-o', pdb_out]
    with tempfolder():
        subprocess.check_call(command)
        with open(pdb_out) as fh:
            assert fh.read().startswith('@<TRIPOS>MOLECULE')


@unittest.skipUnless(_has_program('tleap'), 'Must have tleap')
def test_mutation():
    # mol2
    pdb_fn = get_fn('2igd/2igd.pdb')
    pdb_out = 'out.pdb'

    # no whitespace
    command = [
        'pdb4amber', '-i', pdb_fn, '-o', pdb_out, '-m', '1-ALA,2-ALA,3-ALA'
    ]

    with tempfolder():
        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        assert set(res.name for res in parm.residues[:3]) == {"ALA"}

    # with whitespace
    command = [
        'pdb4amber', '-i', pdb_fn, '-o', pdb_out, '-m', '1-ALA, 2-ALA, 3-ALA'
    ]

    with tempfolder():
        subprocess.check_call(command)
        parm = pmd.load_file(pdb_out)
        assert set(res.name for res in parm.residues[:3]) == {"ALA"}


def test_keep_altlocs():
    ''' e.g: pdb4amber 2igd.pdb --keep-altlocs --reduce'''
    pdb_fn = get_fn('2igd/2igd.pdb')
    command = ['pdb4amber', pdb_fn, '--keep-altlocs', '--reduce']

    with tempfolder():
        output = subprocess.check_output(command).decode()
        input_pdb = StringIO(output)
        input_pdb.seek(0)
        parm = pmd.read_PDB(input_pdb)
        res4 = parm.residues[4]
        for atom in res4.atoms:
            if atom.name.startswith('CB') or atom.name.startswith('CG'):
                assert atom.other_locations


def test_find_gaps():
    expected_line = 'gap of 4.134579 A between MET 1 and PRO 2'
    pdb_fh = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fh)
    parm_gap = parm[':1,3']
    with tempfolder():
        fn = 'gap.pdb'
        parm_gap.save(fn)
        # remove TER to have gap
        subprocess.check_call(
            'cat {} | sed "/TER/d" > new1.pdb'.format(fn), shell=True)
        output = subprocess.check_output(
            ['pdb4amber', 'new1.pdb', '--logfile=stdout']).decode()
        assert expected_line in output


def get_num_ters(fn):
    with open(fn) as fh:
        num_ters = 0
        for line in fh:
            if line.startswith("TER"):
                num_ters += 1
    return num_ters


def test_noter():
    pdb_fh = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fh)
    parm_gap = parm[':1,3']

    with tempfolder():
        fn = 'gap.pdb'
        noter_fn = 'noter.pdb'
        parm_gap.save(fn)
        assert get_num_ters(fn) == 2
        output = subprocess.check_output(
            ['pdb4amber', fn, '--noter', '-o', noter_fn])
        assert get_num_ters(noter_fn) == 0


def test_noter():
    pdb_fh = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fh)
    parm_gap = parm[':1,3']

    with tempfolder():
        fn = 'gap.pdb'
        noter_fn = 'noter.pdb'
        parm_gap.save(fn)
        assert get_num_ters(fn) == 2
        output = subprocess.check_output(
            ['pdb4amber', fn, '--noter', '-o', noter_fn])
        assert get_num_ters(noter_fn) == 0


def test_noter_using_run_method():
    pdb_fh = get_fn('2igd/2igd.pdb')
    parm = pmd.load_file(pdb_fh)
    parm_gap = parm[':1,3']

    with tempfolder():
        fn = 'gap.pdb'
        noter_fn = 'noter.pdb'
        parm_gap.save(fn)
        assert get_num_ters(fn) == 2
        pdb4amber.run(
            noter_fn,
            fn,
            arg_elbow=True,
            arg_reduce=False,
            arg_logfile=sys.stderr,
            arg_conect=False,
            arg_noter=True,
        )
        assert get_num_ters(noter_fn) == 0


def test_increase_code_coverage_for_small_stuff():
    with pytest.raises(RuntimeError):
        pdb4amber.run('fake.pdb', 'fake.pdb')
