import os

files = """avg.pdb
cjt1.dat
dummy.dat
fluct.agr
ired.vec
matrix_ds2.dat
na.pk
noe
nonortho.dx
orderparam
plateau.norm.dat
project.dat
spam.info
test.dat
test.pdb
v0
v0.cjt
v0.cmt
""".split()

for fn in files:
    try:
        os.remove(fn)
    except OSError:
        pass
