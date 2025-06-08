#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Note: To use the 'upload' functionality of this file, you must:
#   $ pip install twine

import io
import os
import sys
from shutil import rmtree

from setuptools import find_packages, setup, Command
from setuptools.command.install import install

# Package meta-data.
NAME = 'packmol_memgen'
DESCRIPTION = 'A PACKMOL wrapper for packing proteins into membranes for AMBER molecular dynamics simulations.'
URL = 'http://cpclab.uni-duesseldorf.de/'
EMAIL = 'schottve@hhu.de'
AUTHOR = 'Stephan Schott-Verdugo'
REQUIRES_PYTHON = '>=3'
VERSION = None

# What packages are required for this module to be executed?
REQUIRED = ['matplotlib', 'pandas']

here = os.path.abspath(os.path.dirname(__file__))

with io.open(os.path.join(here, 'README.rst'), encoding='utf-8') as f:
    long_description = '\n' + f.read()

about = {}
if not VERSION:
    with open(os.path.join(here, NAME, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION

# read arguments to setup script
inside_cmake = False
memembed = True

args_to_pass_on = []

for i, arg in enumerate(sys.argv):
    
    # provided by Amber's CMake build system because it builds packmol and memembed separately
    if arg == '--inside-cmake':
        print('Inside Amber CMake build system, skipping build of native executables')
        inside_cmake = True
    elif arg == '--no-memembed':
        print('Installation will continue without including Memembed. Protein orientation will not be available!')
        memembed = False
    else:
        args_to_pass_on.append(arg)

# replace argument list with our sanitized one
sys.argv = args_to_pass_on

class UploadCommand(Command):
    """Support setup.py upload."""

    description = 'Build and publish the package.'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds…')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building Source and Wheel (universal) distribution…')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPi via Twine…')
        os.system('twine upload dist/*')

        self.status('Pushing git tags…')
        os.system('git tag v{0}'.format(about['__version__']))
        os.system('git push --tags')
        
        sys.exit()

class MakefileInstall(install):
    """Customized setuptools install command for Makefile."""
    def run(self):
        if not inside_cmake:
            if memembed:
                os.system("cd packmol_memgen && make clean; make all")
            else:
                os.system("cd packmol_memgen && make clean; make install.packmol")
        install.run(self)

setup(
    name=NAME,
    version=about['__version__'],
    description=DESCRIPTION,
    long_description=long_description,
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(exclude=('tests',)),
    install_requires=REQUIRED,
    include_package_data=True,
    license='GPL2',
    scripts=['packmol-memgen'],
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Operating System :: Unix',
    ],
    cmdclass={
        'install': MakefileInstall,
        'upload': UploadCommand,
    },
)
