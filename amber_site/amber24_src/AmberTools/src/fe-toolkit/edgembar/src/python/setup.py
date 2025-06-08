#!/usr/bin/env python3

from setuptools import find_packages, setup
from glob import glob

install_requires = [
    'numpy',
    #'matplotlib>=3.5.1'
    'matplotlib'
    ]

package_data = {
    "edgembar": ["pkgdata/*.js",
                 "pkgdata/*.css",
                 "pkgdata/__init__.py"]
    }

scripts = glob("bin/*.py")

setup( name="edgembar",
       version="3.0",
       description="Companion library to the edgembar C++ program " + \
       "for analyzing free energy molecular dynamics simulations",
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
       packages=["edgembar","edgembar.pkgdata"],
       package_dir={"edgembar": "./lib/edgembar"} )

