#!/bin/usr/bash
for dir in ./*
  do
      #ls $dir | wc -l
      fNum="$(ls $dir | wc -l)"
      if [[ "$fNum" -ne 3 ]]; then
          echo $dir
      fi
  done