#!/usr/bin/env python3

from setuptools import find_packages, setup
from glob import glob

install_requires = [
    #'matplotlib>=3.5.1',
    'matplotlib',
    'numpy',
    'scipy'
    ]

# can use 'joblib', but not required

package_data = {}

#scripts = []
scripts = glob("bin/*.py")

setup( name="ndfes",
       version="3.0",
       description="Companion library to the ndfes C++ program " + \
       "for analyzing umbrella window simulations",
       author="Timothy J. Giese",
       author_email="TimothyJGiese@gmail.com",
       platforms=["any"],
       license="MIT",
       url=None,
       python_requires='>3.5',
       install_requires=install_requires,
       include_package_data=True,
       package_data=package_data,
       scripts=scripts,
       packages=["ndfes","ndfes.amber",
                 "ndfes.constants","ndfes.normalmode",
                 "ndfes.gaussian","ndfes.deprecated"],
       package_dir={"ndfes": "./lib/ndfes"} )

