
cannot_run() {
   echo "Warning: $1 cannot be executed !"
   echo "         It is likely that it was not installed correctly."
   echo "         For the parallel case check mpi4py."
   echo "Skipping test."
   exit 0
}
