#!/bin/bash



NUMSRC=$(./NUM_SRC.sh)


for i in `seq -w 0001 "${NUMSRC}"`
do
        ln -s -f -T "../../Par_file_${1}" "run${i}/DATA/Par_file"
done
