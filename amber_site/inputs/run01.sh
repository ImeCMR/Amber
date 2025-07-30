#!/bin/bash
#SBATCH --job-name=eq.job
##SBATCH --output=ala_d_sol.pdb
#SBATCH --error=eq.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=bl.perera@ufl.edu
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --account=alberto.perezant
#SBATCH --qos=alberto.perezant-b
#SBATCH --ntasks=1
date;hostname;pwd

# Set up the environment for your local Amber installation
export AMBERHOME=/orange/alberto.perezant/bl.perera/amber24
export PATH=$AMBERHOME/bin:$PATH
export LD_LIBRARY_PATH=$AMBERHOME/lib:$LD_LIBRARY_PATH

#pdb4amber -i tau_alpha_m01.pdb -o tau_alpha_m01_clean.pdb --reduce --add-missing-atoms -l tau.log
$AMBERHOME/bin/sander -O -i run01.in -o test00.out -p ../setup/test00.prmtop -c ../eq/test00_eq.rst

$AMBERHOME/bin/pmemd -O -i run01.in -o test00.out -p ../setup/test00.prmtop -c ../eq/test00_eq.rst