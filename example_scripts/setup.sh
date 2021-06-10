#!/bin/bash

#this script sets up the directories and links the executables

SPECFEM_DIR="/home/user/specfem2d/"
N_Source=$(./NUM_SRC.sh)

mkdir -p ./DATA
mkdir -p ./OUTPUT_FILES
mkdir -p ./EXECUTABLES

# copy specfem Par_file files
cp "${SPECFEM_DIR}/DATA/Par_file" ./DATA
cp "${SPECFEM_DIR}/DATA/SOURCE" ./DATA
cp "${SPECFEM_DIR}/EXAMPLES/semi_infinite_homogeneous/DATA/interfaces_elastic_analytic.dat" ./DATA/interfaces.dat



rm -fr ./EXECUTABLES/xmeshfem2D
rm -fr ./EXECUTABLES/xspecfem2D
ln -s "${SPECFEM_DIR}/bin/xmeshfem2D" ./EXECUTABLES/
ln -s "${SPECFEM_DIR}/bin/xspecfem2D" ./EXECUTABLES/


for i in `seq -w 0001 $N_Source`
do
  mkdir -p "run${i}"
  mkdir -p "run${i}/DATA"
  mkdir -p "run${i}/EXECUTABLES"
  mkdir -p "run${i}/OUTPUT_FILES"
  mkdir -p "run${i}/SEM"

  rm -fr "run${i}/EXECUTABLES/xmeshfem2D"
  rm -fr "run${i}/EXECUTABLES/xspecfem2D"
  ln -s "${SPECFEM_DIR}/bin/xmeshfem2D" "run${i}/EXECUTABLES"
  ln -s "${SPECFEM_DIR}/bin/xspecfem2D" "run${i}/EXECUTABLES"

done
