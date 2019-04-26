#!/bin/usr/bash

### Check if any of the outputs do not have adequate number of files

for dir in ./*
  do
      #ls $dir | wc -l
      fNum="$(ls $dir | wc -l)"
      if [[ "$fNum" -ne 3 ]]; then
          echo $dir
      fi
  done