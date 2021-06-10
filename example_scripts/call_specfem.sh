#!/bin/bash

## modify mpi options as necessary


#run forward simulation in specfem2d
echo
echo `date`
echo
N_Source=$(./NUM_SRC.sh)


# Get the number of processors
NPROC=`grep ^NPROC DATA/Par_file | cut -d = -f 2 | cut -d \# -f 1 | tr -d ' '`

echo
echo "Number of processors: $NPROC"
echo

if [ "$NPROC" -eq 1 ]; then
  #this is a serial simulation
  echo
  echo "running mesher"
  echo
  ./EXECUTABLES/xmeshfem2D
else
  #this is an MPI simulation
  echo "runnning mesher on $NPROC processors..."
  echo
  mpirun -np $NPROC --map-by node --bind-to none --mca mpi_cuda_support 0 --mca mca_base_component_show_load_errors 0 ./EXECUTABLES/xmeshfem2D > /dev/null
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

# runs simulation
if [ "$NPROC" -eq 1 ]; then
  # This is a serial simulation
  echo
  echo "running solver..."
  echo
  ./EXECUTABLES/xspecfem2D
else
  # This is a MPI simulation
  echo
  echo "running solver on $NPROC processors..."
  echo
  mpirun -np $NPROC --map-by node --bind-to none --mca mpi_cuda_support 0 --mca mca_base_component_show_load_errors 0 ./EXECUTABLES/xspecfem2D > /dev/null
fi
# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi


echo
echo "done"
echo `date`
