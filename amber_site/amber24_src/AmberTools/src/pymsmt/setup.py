
from setuptools import setup

# Packages in python MSMT toolbox
packages = ['pymsmt', 'pymsmt.api', 'pymsmt.mcpb', 'pymsmt.mol', 'pymsmt.ipmach']

# Modules
modules = ['pymsmt.exp', 'pymsmt.title', 'pymsmt.lib']

# Scripts
scripts = ['pymsmt/tools/MCPB.py', 'pymsmt/tools/OptC4.py', 'pymsmt/tools/PdbSearcher.py',
           'pymsmt/tools/espgen.py', 'pymsmt/tools/CartHess2FC.py', 'pymsmt/tools/IPMach.py',
           'pymsmt/tools/car_to_files.py', 'pymsmt/tools/ProScrs.py', 'pymsmt/tools/mol2rtf.py',
           'pymsmt/tools/amb2chm_psf_crd.py', 'pymsmt/tools/amb2chm_par.py',
           'pymsmt/tools/amb2gro_top_gro.py', 'pymsmt/tools/metalpdb2mol2.py']

if __name__ == '__main__':

    setup(name='pyMSMT',
          version='22.0', # For AmberTools 22
          description='python Metal Site Modeling Toolbox',
          author='Pengfei Li',
          author_email='ldsoar1990 -at- gmail.com',
          license='GPL v3 or later',
          packages=packages,
          py_modules=modules,
          scripts=scripts)
