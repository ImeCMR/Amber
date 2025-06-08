#!/bin/sh

cat > __fix__.awk << EOF
# Get rid of the automake configure dependencies
/^am__configure_deps/ {getline; next}

# Get rid of configure 
/^\\\$\\(top_srcdir\\)\\/configure/ {getline; next}

# Remove Makefile.am from the DIST_COMMON distfile list
/\\\$\\(srcdir\\)\\/Makefile\\.am/ {
   myline = \$0
   sub(/\\\$\\(srcdir\\)\\/Makefile\\.am/, "", myline)
   print myline
   next
}

# Remove dependencies of Makefile.in
/^\\\$\\(srcdir\\)\\/Makefile\\.in/ {
#  print "FOUND MY LINE!"
   getline
   while (match(\$0, /^[\\t\\r\\n]/)) {getline}
   print \$0
   next
}

/.*/ {print \$0}
EOF

for d in netcdf-4.3.0 netcdf-fortran-4.4.4; do
   for f in `find $d -name "Makefile.in"`; do
   
      echo "Fixing $f"
   
      awk -f __fix__.awk $f > tmp; /bin/mv tmp $f
   
   done
done

/bin/rm -f __fix__.awk
