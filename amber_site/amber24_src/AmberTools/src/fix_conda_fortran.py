# This script fixes libgfortran issues in conda if they exist
import os
import sys
import subprocess

if 'linux' not in sys.platform:
    # Not sure what to do on other OSes, so don't do anything
    sys.exit(0)

if os.getenv('AMBERHOME') is None:
    sys.stderr.write('AMBERHOME must be set to try to fix libgfortran issues\n')
    sys.exit(0)

# Go to $AMBERHOME, since that's where it's easiest to work from
if not os.path.isdir(os.getenv('AMBERHOME')):
    sys.stderr.write('AMBERHOME [%s] is not a directory!\n' %
                     os.getenv('AMBERHOME'))
os.chdir(os.getenv('AMBERHOME'))

# This script only fixes Amber's conda
pyexe = os.path.realpath(os.path.abspath(sys.executable))
ambhome = os.path.realpath(os.path.abspath(os.getenv('AMBERHOME')))

if not pyexe.startswith(ambhome):
    sys.stderr.write("I only work on Amber's miniconda\n")
    sys.exit(0)

# Now we know we are using Amber's miniconda... we should be able to import
# pytraj. If not, check that the error message says something about libgfortran
# and GFORTRAN, which is the well-known error reported here:
# https://github.com/ContinuumIO/anaconda-issues/issues/686. Once Continuum IO
# fixes this, this script will do nothing (since importing pytraj will be
# successful)

try:
    import pytraj
except ImportError as e:
    if not 'libgfortran' in str(e) and not 'GFORTRAN' in str(e):
        # This isn't the bug we are here to catch
        sys.exit(0)
    # Find what version of libgfortran.so.x.x.x we are looking for
    import re
    try:
        libgfortran_version = re.findall(r'libgfortran.so[\.0-9]+', str(e))[0]
    except IndexError:
        # We couldn't find the libgfortran version...
        sys.stderr.write('Could not find libgfortran version to replace...\n')
        sys.exit(0)
else:
    # Nothing broken!!
    sys.exit(0)

print('Attempting to fix libgfortran breakage from Miniconda...')

# Assume the compiler is simply gfortran
for line in subprocess.check_output(['gfortran', '-print-search-dirs']).decode().split('\n'):
    if not line.startswith('libraries'):
        continue
    for folder in line[line.index('=')+1:].split(':'):
        testname = os.path.join(folder, libgfortran_version)
        if os.path.exists(testname):
            system_file = os.path.realpath(os.path.abspath(testname))
            break
    else:
        sys.stderr.write('Could not find %s in gfortran standard search paths!\n'
                         % libgfortran_version)
        sys.exit(0)
    break # Found our system_file
else:
    sys.stderr.write('Could not find standard search directories for libraries\n')
    sys.exit(0)

# OK... we have our system compiler libgfortran.so file now. Go ahead and drop
# it in place
fname = os.path.split(system_file)[1] # Get the base name
import shutil
# Back up the existing libgfortran
shutil.copy(
        os.path.join('miniconda', 'lib', fname),
        os.path.join('miniconda', 'lib', fname+'.bak')
)
shutil.copy(system_file, os.path.join('miniconda', 'lib', fname))
