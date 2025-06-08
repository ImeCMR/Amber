#!/usr/bin/env python

from setuptools import setup
import os
import shutil

# First the ParmedTools packages:
packages = ['MMPBSA_mods']

# Scripts
scripts = ['MMPBSA.py', 'ante-MMPBSA.py', 'MMPBSA.py.MPI']

if __name__ == '__main__':

    shutil.copyfile('MMPBSA.py', 'MMPBSA.py.MPI')
    shutil.copymode('MMPBSA.py', 'MMPBSA.py.MPI')
    setup(name='MMPBSA.py',
          version='16.0',
          description='Program for carrying out MM/PBSA-like end-state '
                      'free energy calculations',
          author='Jason Swails, T. Dwight McGee Jr., and Bill R. Miller III',
          author_email='jason.swails -at- gmail.com',
          url='http://ambermd.org',
          license='GPL v2 or later',
          packages=packages,
          scripts=scripts)

    # Delete our temporary file
    os.unlink('MMPBSA.py.MPI')
