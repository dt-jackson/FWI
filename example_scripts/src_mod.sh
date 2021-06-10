#!/bin/bash


f_n_old=".\/DATA\/src_file.txt"
initial_loc=-999

N_Source=$(./NUM_SRC.sh)

for i in `seq -w 0001 $N_Source`
do
  loc=$(bc<<<"0+(1*($i-1))")
  f_n_new=".\/run"$i"\/DATA\/src_file.txt"
  #echo $f_n_new
  #echo $loc
  sed -i "s/${initial_loc}/${loc}/g" "./run${i}/DATA/SOURCE"
  sed -i "s/${f_n_old}/${f_n_new}/g" "./run${i}/DATA/SOURCE"
done
