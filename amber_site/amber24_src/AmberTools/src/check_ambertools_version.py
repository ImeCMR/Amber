import sys
import os

# change to AMBERHOME
os.chdir(os.path.dirname(__file__) + '/../../')

# files : List[Tuple[filename, List[expected_string]]]
file_combo = [
    ('AmberTools/src/sander/multisander.F90', [
        '!        SANDER, version {}',
        'character (len=4), parameter :: VERSION_STRING = "{}.0"',
        ]),
    ('AmberTools/src/configure2', [
        'echo "       AMBER {} requires CUDA version 7.5 or 8.0"',
        ]),
    ('AmberTools/src/etc/setup.py', [
        "version='{}.0',",
        ]),
]

expected_version = sys.argv[1]

n_errors = 0

for fn, expected_lines in file_combo:
    print("checking {}".format(fn))
    with open(fn) as fh:
        LINES = [line.strip() for line in fh.readlines()]
        for line in expected_lines:
            expected_string = line.format(expected_version)
            if not expected_string in LINES:
                n_errors += 1
                print("Not found: '{}' in {}".format(expected_string, fn))
                

if n_errors == 0:
    print("Sounds good")
else:
    print("WARNING you")
    sys.exit(1)
