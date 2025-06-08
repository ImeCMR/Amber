#!/bin/sh

if perl < /dev/null > /dev/null 2>&1  ; then
    export PERL5LIB=$AMBER_SOURCE/AmberTools/src/FEW/additional_libs/PerlMol-0.3500
    cd $AMBER_SOURCE/AmberTools/src/FEW/additional_libs/PerlMol-0.3500
    perl Makefile.PL PREFIX="$AMBER_PREFIX/lib/perl/"
    make
    make install
    if [ ! -d $AMBER_PREFIX/lib/perl/FEW_libs ]; then
        mkdir "$AMBER_PREFIX/lib/perl/FEW_libs"
    fi
    cp $AMBER_SOURCE/AmberTools/src/FEW/libs/*pm $AMBER_PREFIX/lib/perl/FEW_libs/
    cp $AMBER_SOURCE/AmberTools/src/FEW/additional_libs/*pm $AMBER_PREFIX/lib/perl/FEW_libs/
    cp -rf $AMBER_SOURCE/AmberTools/src/FEW/examples/command_files $AMBER_PREFIX/lib/perl/FEW_libs/
    cp -rf $AMBER_SOURCE/AmberTools/src/FEW/examples/input_info $AMBER_PREFIX/lib/perl/FEW_libs/
    cd $AMBER_SOURCE/AmberTools/src/FEW
else
	echo "Error:  No Perl available on system. Please install Perl !"
	exit 1
fi
