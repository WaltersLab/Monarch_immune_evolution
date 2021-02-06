#!/bin/bash

#### This is a script for running ANGSD to get popgen stats with several inputs ####

## Load angsd 
module load angsd/0.919

## Info of command line args
# $1 = population name for bam.filelist 
# $2 = gene name (ID)
# $3 = "region" 
# $4 = sites (0f or 4f) 

## Set directories 
dir=''
dir_out=''
dir_r=''

## Running ANGSD to calculate popgen stats 
# SFS pre-generated 

## step 4 Calculate the thetas
angsd sites index $dir/Sites/$2"_sites_"$4.txt        # index the individual gene site file 
angsd -bam $dir/$1#_bam.filelist -sites $dir/Sites/$2"_sites_"$4.txt \
-r $3 -out $2"_"$4"_outFold" -doThetas 1 -doSaf 1 -pest $1"_"$4"_outFold".sfs \
-anc $dir_r/Dp_genome_v3_masked.fasta -GL 1 -fold 1

## step 5 Estimate Tajimas D
thetaStat do_stat $2"_"$4"_outFold".thetas.idx 

## step 6 Print out ThetaStat values 
thetaStat print $2"_"$4"_outFold".thetas.idx > $2"_"$4"_outFold".thetas.idx.txt

## sort outputs
cp $2"_"$4"_outFold.thetas.idx.txt" $dir_out/All_stats/$1  # copy out the pi and theta data
cp $2"_"$4"_outFold.thetas.idx.pestPG" $dir_out/All_stats/$1   # copy out the Tjd data
mkdir $dir_out/$1/out_$2
mv $2* $dir_out/$1/out_$2/  # sort all files into their own folder 
