#!/bin/usr/bash

### Check if any genes do not have adequate number of output files ### 

for dir in ./*
  do
      #ls $dir | wc -l
      fNum="$(ls $dir | wc -l)"
      if [[ "$fNum" -ne 5 ]]; then
          echo $dir
      fi
  done