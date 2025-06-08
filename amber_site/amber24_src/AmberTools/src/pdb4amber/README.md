Usage
-----
```bash
$ pdb4amber --help

# Some examples
# simplest case, print output pdb
pdb4amber my.pdb 

# simplest case, save output to file
pdb4amber my.pdb -o out.pdb

# compat mode
pdb4amber -i my.pdb -o out.pdb

# pipe
cat my.pdb | pdb4amber -o out.pdb

# save to different formats (e.g: mol2)
pdb4amber my.pdb -o out.mol2

# use `reduce` program to add hydgron
pdb4amber my.pdb --reduce -o out.pdb

# process other formats (e.g: .cif)
pdb4amber my.cif --reduce -o out.pdb

# process from URL
pdb4amber https://raw.githubusercontent.com/ParmEd/ParmEd/master/test/files/4LZT.cif -o out.pdb

# fetch structure by its pdbid and process
pdb4amber 1tsu --pdbid --reduce -o out.pdb

# logfile
pdb4amber my.pdb -o out.pdb --logfile=my.log
pdb4amber my.pdb -o out.pdb --logfile=stdout
```

