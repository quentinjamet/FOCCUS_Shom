#!/bin/bash

for iyr in {2019..2019}; do
  for imonth in {02..04}; do
    for imem in {000..016}; do
      echo "year: " ${iyr}  "; month: " ${imonth} "; member: " ${imem}
      fileList=( $(ls ./${iyr}/${imonth}/mem${imem}/*oce*.nc) )
      for fileIN in ${fileList[@]}; do
        fileOUT=$( echo ${fileIN} | sed -e 's/\oce/\zos/g')
        ncks -v zos ${fileIN} ${fileOUT}
      done
    done
  done
done
