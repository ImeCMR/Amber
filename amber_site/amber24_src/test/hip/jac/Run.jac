#!/bin/csh -f
#TEST-PROGRAM pmemd.hip
#TEST-DESCRIP TO_BE_DEtermined
#TEST-PURPOSE regression, basic
#TEST-STATE   undocumented

#$1 = PREC_MODEL
#$2 = NETCDF

if( ! $?DO_PARALLEL ) then
  setenv DO_PARALLEL " "
  if( $?TESTsander ) then
      set sander = $TESTsander
  else
      set sander = ${AMBERHOME}/bin/pmemd.hip_$1
  endif
else
  if( $?TESTsander ) then
      set sander = $TESTsander
  else
      set sander = ${AMBERHOME}/bin/pmemd.hip_$1.MPI
  endif
endif

cat > mdin <<EOF
 short md, jac, power 2 FFT
 &cntrl
   ntx=7, irest=1,
   ntc=2, ntf=2, tol=0.0000001, 
   nstlim=20, 
   ntpr=5, ntwr=10, 
   dt=0.001, 
   cut=9.,
   ntt=0, temp0=300., ig=71277,
 /
 &ewald
  nfft1=64,nfft2=64,nfft3=64,netfrc=0,
 /
EOF

set output = jac.out

set output_save = $output.GPU_$1

touch dummy
$DO_PARALLEL $sander -O -o $output -c inpcrd.equil -r restrt -x mdcrd < dummy || goto error

#Use different dacdif tolerances based on precision model.
if ( "$1" == "DPFP" ) then
  #7 sig figs
  ../../dacdif -r 1.0e-7 $output_save $output
else
  #3 sig figs
  ../../dacdif -r 3.1e-3 $output_save $output
endif

/bin/rm -f mdin logfile mdinfo dummy restrt mdcrd
exit(0)

error:
echo "  ${0}:  Program error"
exit(1)










