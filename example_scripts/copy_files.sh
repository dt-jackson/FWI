#!/bin/bash

# The Par_file should be modified so it is correct before running this script


# copy parameter files from DATA directory 
cp ./DATA/Par_file ./Par_file_f
sed -e 's/SIMULATION_TYPE *\t* = *\t* 3/SIMULATION_TYPE = 1/g' DATA/Par_file > tmp_par_file
sed -e 's/SAVE_FORWARD *\t* = *\t* .false./SAVE_FORWARD = .true./g' tmp_par_file > Par_file_f

sed -e 's/SIMULATION_TYPE *\t* = *\t* 1/SIMULATION_TYPE = 3/g' Par_file_f > tmp_par_file
sed -e 's/SAVE_FORWARD *\t* = *\t* .true./SAVE_FORWARD = .false./g' tmp_par_file > Par_file_a
rm -f tmp_par_file


#copy source and interface files from DATA directory
N_Source=$(./NUM_SRC.sh)

for i in `seq -w 0001 $N_Source`
do
  cp -p ./DATA/SOURCE "./run${i}/DATA"
  cp -p ./DATA/interfaces.dat "./run${i}/DATA"

  ### uncomment this line and add correct file name 
  ### to copy a source-time function file as well
  #cp -p ./DATA/src_file.txt "./run${i}/DATA" 
done
