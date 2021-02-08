import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
import re
import random 
from Bio import SeqIO
from Bio.Seq import Seq
import argparse

 

 
 
def get_args():
    parser = argparse.ArgumentParser(description='''This script gives nucleotide diversity , watterson's theta  and Tajima's D  from allele frequency files. 
        For this you will need to have the follwoing files in a directory
        1) All positions for each codon position (codon1: codon_df_1.csv.gz; codon2: codon_df_2.csv.gz; codon3: codon_df_3.csv.gz; 4fold: codon_df_4d.csv.gz; 0fold: codon_df_0d.csv.gz;)
        2) for intergenic and all sites you will need to provide the refrence genome file as refrence.fa
        -F allele frequncy file (from VCF tools)
        -S input sites (CSV file with sites to be used)
        -I sites type (codon1/codon2/codon3/4fold/0fold)
        -P Path
        -Sec Section
        ''')

    parser.add_argument('-F', '--Freq', type=str, required=True)
    parser.add_argument('-S', '--input_sites', type=str, required=True)
    parser.add_argument('-I', '--sites_type', type=str, required=True)
    parser.add_argument('-P', '--Path', type=str, required=True)
    parser.add_argument('-Sec', '--part', type=int, required=True)
    return parser.parse_args()

args = get_args()


#path='/scratch/vt20265/bam_files'
#frequency_file='/scratch/vt20265/1977.frq'
#input_sites='codon_df_1.csv'
#sites_type='codon1'
#part=0



pwd=args.Path
os.chdir(pwd)

frequency_file=args.Freq
input_sites=args.input_sites
sites_type=args.sites_type
part=int(args.part)
condons=['4fold','0fold','codon1','codon2','codon3']



chromosome=pandas.read_csv("/scratch/vt20265/bam_files/windows.bed", sep='\t', names=['CHROM','BIN_START','BIN_END'], header=None)
temp=np.array_split(chromosome, 20)[part]
temp=list(set(temp.CHROM))
chromosome=chromosome[chromosome['CHROM'].isin(temp)]



def join_data_bases(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END','POS'], how='outer')
        nextone=new_2
    return nextone


def fasta_dict(fasta):
    fasta = SeqIO.parse(fasta,"fasta")
    seq_dict = {}
    for record in fasta:
        seq_dict[record.id]=record.seq
    return seq_dict



def afs_basic_stats(afs):
    ''' 
    Calculate Tajima's D,theta watterson and nucleotide diversity from the allele frequency spectrum.
    Input: the site frequency spectrum as a list for x chromosmes:
    [1/x, 2/x, 3/x....x-1/x, total valid sites]
    It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list. 
    It returns, pi, theta watterson and TajD
    Example:
    >>> afs_basic_stats([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[1]: (0.093333333333333338, 0.08759124087591241, 0.0020911421903314943)
    >>> afs_basic_stats([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[2]: (0.093333333333333338, 0.08759124087591241, 0.0020911421903314943)
    '''
    Sites = afs_thetaW(afs)[1]
    mythetaW = afs_thetaW(afs)[0] # calculate Theta Watterson
    mypi = afs_nucl_diversity(afs) # calculate nucleotide diversity
    varD =  variance_D(afs) # calculate the variance for d
    TajD = ((mypi*Sites)- (mythetaW*Sites)) / varD # calculate Tajima's D
    print (" returning pi, theta and TajD....")
    return mypi, mythetaW, TajD


def afs_thetaW(afs):
    ''' 
    Calculate theta watterson from the allele frequency spectrum.
    Input: the site frequency spectrum as a list for x chromosmes:
    [1/x, 2/x, 3/x....x-1/x, total valid sites]
    It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list.
    Example:
    >>> afs_thetaW([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[1]: 0.08759124087591241
    >>> afs_thetaW([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[2]: 0.08759124087591241
    '''
    nchrom = len(afs) # the number of chromomosomes
    seg_site = sum([int(x) for x in afs[:-1]]) # number of segregating sites
    n_total = float(afs[-1]) # total number of sites
    harmonic_mean = sum([1/float(x) for x in range(1,nchrom)]) # harmonic mean 
    theta_watterson =[(seg_site/harmonic_mean) / n_total, n_total]
    print ("assuming",nchrom,"chromosomes",seg_site,"variable sites")
    return theta_watterson


def afs_nucl_diversity(afs):
    ''' 
    Calculate nucleotide diveristy (pi) from the allele frequency spectrum.
    Input: the site frequency spectrum as a list for x chromosmes:
    [1/x, 2/x, 3/x....x-1/x, total valid sites]
    It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list.
    Example:
    >>> afs_nucl_diversity([10,9,8,7,6,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[1]: 0.093333333333333338
    >>> afs_nucl_diversity([16,16,8,0,0,200]) ## unfolded spectrum for si chromosomes with 200 valid sites
    Out[2]: 0.093333333333333338
    '''
    n = len(afs) # number of chromosomes
    sf = np.array([float(x) for x in afs[:-1]]) # allelc frequency spectrum
    tot = afs[-1] # number of valid sites
    pi = sf*2* np.array(range(1,n)) * np.array(range(n-1,0,-1)) 
    pi2 = (sum(pi /(n*(n-1)))) /float(tot)
    return pi2


def variance_D(afs):
    '''Calculate Tajima's D,theta watterson and nucleotide diversity from the allele frequency spectrum.
    Input: the site frequency spectrum as a list for x chromosmes:
    [1/x, 2/x, 3/x....x-1/x, total valid sites]
    It works for both folded and unfolded site frequency spectrum but has to define x to -1 chromosomes in the list. 
    '''
    n = len(afs) # number of chromosomes
    a1 = sum([1.0/x for x in range(1,n)])
    b1 = (n+1.)/(3*(n-1)) 
    c1 = b1-( 1/a1)
    e1 = c1/a1

    a2 = sum([1.0/(x**2.) for x in range(1,n)])
    b2 = ( 2*(n**2.+n+3.)  )/ ( (9.0*n) *  (n-1.)      )
    c2 = b2 - ((n+2.)/(a1*n)) + (a2)/(a1**2.)
    e2 = c2 /((a1**2.)+(a2))
    
    s = sum([int(x) for x in afs[:-1]]) # number of segregating sites
    varD = (((e1*s)+(e2*s*(s-1)))**0.5)
    return varD



def get_site_freq_spectrum(goods,codon_df_1):
    codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    temp=pandas.DataFrame({'CHROM':['Fake']*8, 'Ind':range(1,9),'BIN_START':range(1,9),'POS':range(1,9),'major_allele':range(1,9),'major_freq':range(1,9),'minor_allele':range(1,9),'minor_freq':range(1,9),'site_pi':range(1,9),'BIN_END':range(1,9)})
    codon_df_1_temp=pandas.concat([codon_df_1_temp,temp], ignore_index=True, sort=True)
    codon_df_1_temp.minor_freq=codon_df_1_temp.minor_freq.round(2)
    spectrum=codon_df_1_temp.groupby('Ind')['Ind'].count()
    spectrum=spectrum-1
    return spectrum


def get_site_unfolded_freq_spectrum_all(goods,codon_df_1):
    if codon_df_1 == 'all':
        codon_df_1_temp=goods
    else:
        codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    codon_df_1_temp.Derived_freq =codon_df_1_temp.Derived_freq.round(2)
    spectrum=codon_df_1_temp.groupby('Derived_freq')['Derived_freq'].count()
    return spectrum

def get_site_unfolded_freq_spectrum_2(goods,codon_df_1):
    if codon_df_1 != 'all':
        codon_df_1_temp=pandas.merge(goods,codon_df_1, on=['CHROM','POS'], how='inner')
    else:
        codon_df_1_temp=goods
    temp=pandas.DataFrame({'CHROM':['Fake']*29, 'POS':range(0,29), 'Derived_freq':[0.0, 0.04, 0.07, 0.11, 0.14, 0.18, 0.21, 0.25, 0.29, 0.32, 0.36, 0.39, 0.43, 0.46, 0.5, 0.54, 0.57, 0.61, 0.64, 0.68, 0.71, 0.75, 0.79, 0.82, 0.86, 0.89, 0.93, 0.96, 1.0]})
    codon_df_1_temp=pandas.concat([codon_df_1_temp,temp], ignore_index=True, sort=True)
    codon_df_1_temp.Derived_freq =codon_df_1_temp.Derived_freq.round(2)
    spectrum=codon_df_1_temp.groupby('Derived_freq')['Derived_freq'].count()
    return spectrum



def output_good_sites_fequency_new(freq_file1,chr_filter):
    header=['CHROM', 'POS', 'N_ALLELES', 'N_CHR', 'al_1_', 'al_2_','al_3_','al_4_']
    a=pandas.read_table(freq_file1,  names=header, engine='c',skiprows=1, error_bad_lines=False)
    merged=a[(a['N_CHR'] >= int(chr_filter))]
    bads=[]
    goods=[]
    new_output=pandas.DataFrame(columns=['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq'])
    temp2=merged[merged['N_ALLELES']<=2]
    temp2['minor_allele']=''
    temp2['minor_freq']=0
    temp2['minor_allele'][temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float)]=temp2.al_2_.str[0]
    temp2['minor_allele'][temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float)]=temp2.al_1_.str[0]
    temp2['minor_freq'][temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float)]=temp2.al_2_.str[2:].astype(float)
    temp2['minor_freq'][temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float)]=temp2.al_1_.str[2:].astype(float)
    #
    temp2['major_allele']=''
    temp2['major_freq']=0
    #temp2.assign(major_allele=(temp2.al_1_.str[0]).where(float(temp2.al_1_.str[2:)) > float(temp2.al_2_.str[2:))), 0))
    temp2['major_allele'][temp2.al_2_.str[0] == temp2.minor_allele ]=temp2.al_1_.str[0]
    temp2['major_freq'][temp2.al_2_.str[0] == temp2.minor_allele ]=temp2.al_1_.str[2:].astype(float)
    temp2['major_allele'][temp2.al_1_.str[0] == temp2.minor_allele ]=temp2.al_2_.str[0]
    temp2['major_freq'][temp2.al_1_.str[0] == temp2.minor_allele ]=temp2.al_2_.str[2:].astype(float)
    temp2=temp2[['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq']]
    temp2=temp2[temp2.minor_freq>=0]
    new_output=new_output.append(temp2, ignore_index=True)
    temp2=merged[merged['N_ALLELES']==3]
    temp2['minor_allele']=''
    temp2['minor_freq']=0
    temp2['major_allele']=''
    temp2['major_freq']=0
    temp2['minor_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['major_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['minor_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['major_allele'][ (temp2.al_3_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[0]
    #
    temp2['minor_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['minor_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_3_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_2_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    #
    temp2['minor_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    temp2['major_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['minor_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[0]
    temp2['major_allele'][ (temp2.al_2_.str[2:].astype(float) == 0.0) & (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    #
    temp2['minor_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['minor_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_1_.str[2:].astype(float)
    temp2['major_freq'][ (temp2.al_2_.str[2:].astype(float) == 0.0) &   (temp2.al_1_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    #
    temp2['minor_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    temp2['major_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['minor_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[0]
    temp2['major_allele'][ (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[0]
    #
    temp2['minor_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    temp2['major_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) >= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['minor_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_2_.str[2:].astype(float)
    temp2['major_freq'][  (temp2.al_1_.str[2:].astype(float) == 0.0) &   (temp2.al_2_.str[2:].astype(float) <= temp2.al_3_.str[2:].astype(float))]=temp2.al_3_.str[2:].astype(float)
    #
    temp2=temp2[(temp2['minor_allele']!='') & (temp2['major_allele']!='')]
    temp2=temp2[['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq']]
    temp2['POS']=temp2['POS'].astype('int')
    temp2['BIN_START']=0
    temp2['BIN_START']=(np.floor(temp2['POS']/100000)*100000)+1
    temp2=temp2[temp2.minor_freq>=0]
    new_output=new_output.append(temp2, ignore_index=True)
    #temp2=merged[merged['N_ALLELES'] ==4 ]
    #if len(temp2) > 0:
    #   for row in temp2.itertuples():
    #       allele={}
    #       freq_list=[float(row.al_1_[2:]), float(row.al_2_[2:]),float(row.al_3_[2:]), float(row.al_4_[2:])]
    #       if 0.0 in freq_list:
    #           if float(row.al_1_[2:]) not in [0.0]: allele[row.al_1_[0]]=float(row.al_1_[2:])
    #           if float(row.al_2_[2:]) not in [0.0]: allele[row.al_2_[0]]=float(row.al_2_[2:])
    #           if float(row.al_3_[2:]) not in [0.0]: allele[row.al_3_[0]]=float(row.al_3_[2:])
    #           if float(row.al_4_[2:]) not in [0.0]: allele[row.al_4_[0]]=float(row.al_4_[2:])
    #       else:
    #           bads.append(row)
    #       if  ((len(allele.values()) == 2) and min(allele.values()) != 0.0):
    #           temp={}
    #           temp['CHROM']=[row.CHROM]
    #           temp['POS']=[row.POS]
    #           temp['minor_allele']=[min(allele, key=lambda k: allele[k])]
    #           temp['minor_freq']=[allele[min(allele, key=lambda k: allele[k])]]
    #           temp['major_allele']=[max(allele, key=lambda k: allele[k])]
    #           temp['major_freq']=[allele[max(allele, key=lambda k: allele[k])]]
    #           temp_df = pandas.DataFrame(temp)
    #           new_output=new_output.append(temp_df, ignore_index=True)
    #new_output['BIN_START']=(np.floor(new_output['POS']/100000)*100000)+1
    #new_output['BIN_END']=new_output['BIN_START']+(100000-1)
    new_output['site_pi']=0.0
    total_comb=(int(chr_filter)*(int(chr_filter)-1.0))/2.0
    new_output['site_pi']=((new_output.minor_freq*int(chr_filter))*(new_output.major_freq*int(chr_filter)))/total_comb
    new_output['BIN_END']=0
    new_output['BIN_END']=new_output['BIN_START']+(100000-1)
    return new_output


if sites_type=='4fold':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep=' ',header=None)
elif sites_type=='0fold':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep=' ',header=None)
elif sites_type=='codon1':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep=' ',header=None)
elif sites_type=='codon2':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep=' ',header=None)
elif sites_type=='codon3':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep=' ',header=None)
elif sites_type=='introns':
    all_sites=pandas.read_csv(input_sites,names=['CHROM','POS'], sep='\t',header=None)
elif sites_type=='intergenic':
    all_sites_2=pandas.read_csv(input_sites,skiprows=1,names=['CHROM','BIN_START','BIN_END','N_sites'])
    fasta_DPlex=personal_popgen.fasta_dict('/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa')
    codon_df_1=pandas.read_csv('codon_df_1.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_1=codon_df_1[codon_df_1['CHROM'].isin(temp)]
    codon_df_2=pandas.read_csv('codon_df_2.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_2=codon_df_2[codon_df_2['CHROM'].isin(temp)]
    codon_df_3=pandas.read_csv('codon_df_3.csv',names=['CHROM','POS'], sep=' ',header=None)
    codon_df_3=codon_df_3[codon_df_3['CHROM'].isin(temp)]
    introns=pandas.read_csv('introns.csv',names=['CHROM','POS'], sep=' ',header=None)
    introns=introns[introns['CHROM'].isin(temp)]
    intergenic_exclude=pandas.concat([codon_df_1[['CHROM','POS']], codon_df_2[['CHROM','POS']], codon_df_3[['CHROM','POS']], introns[['CHROM','POS']]], ignore_index=True)
    intergenic_exclude=intergenic_exclude[['CHROM','POS']]
    intergenic=pandas.DataFrame(columns=['CHROM','POS'])
    for i in temp:
        chrom_len=len(fasta_DPlex[i])
        all_site=list(range(1,chrom_len+1))
        exlude_sites=list(intergenic_exclude[(intergenic_exclude.CHROM == i )]['POS'])
        temp_1=pandas.DataFrame({'CHROM':[i]*len(all_site), 'POS':all_site})
        temp_1=temp_1[~temp_1['POS'].isin(exlude_sites)]
        intergenic=intergenic.append(temp_1, ignore_index=True)
    all_sites=intergenic
    all_sites['POS']=all_sites['POS'].astype('int')
    all_sites['BIN_START']=0
    all_sites['BIN_START']=(np.floor(all_sites['POS']/10000)*10000)+1
    all_sites['BIN_END']=all_sites['BIN_START']+(10000-1)
if sites_type=='all':
    all_sites_2=pandas.read_csv(input_sites,skiprows=1,names=['CHROM','BIN_START','BIN_END','N_sites'])


if sites_type in ['4fold','0fold','codon1','codon2','codon3','introns']:
    all_sites['POS']=all_sites['POS'].astype('int')
    all_sites['BIN_START']=0
    all_sites['BIN_START']=(np.floor(all_sites['POS']/10000)*10000)+1
    all_sites['BIN_END']=all_sites['BIN_START']+(10000-1)
    all_sites_2=all_sites.groupby(['CHROM','BIN_START','BIN_END'])['POS'].count()
    all_sites_2=all_sites_2.reset_index()
    all_sites_2=all_sites_2.rename(index=str, columns={'POS': 'N_sites' })


East_freq = output_good_sites_fequency_new(frequency_file, 16)
East_freq = East_freq[East_freq['CHROM'].isin(temp)]
East_freq = East_freq[East_freq['major_allele'].isin(['A','T','G','C'])]
East_freq = East_freq[East_freq['minor_allele'].isin(['A','T','G','C'])]
East_freq['POS']=East_freq['POS'].astype('int')
East_freq['BIN_START']=0
East_freq['BIN_START']=(np.floor(East_freq['POS']/10000)*10000)+1
East_freq['BIN_END']=East_freq['BIN_START']+(10000-1)
East_freq['Ind']=East_freq.minor_freq*16
East_freq['Ind']=East_freq['Ind'].round()
East_freq=East_freq[East_freq.Ind>=1.0]
if sites_type != 'all':
    East_freq=East_freq.merge(all_sites, on=['CHROM','POS','BIN_START', 'BIN_END'], how='inner')


windows=chromosome
East_freq_1=East_freq[East_freq.Ind==1.0]
East_freq_1=East_freq_1.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_1=East_freq_1.reset_index()
East_freq_1=East_freq_1.rename(index=str, columns={'Ind': 'IND_1' })
East_freq_2=East_freq[East_freq.Ind==2.0]
East_freq_2=East_freq_2.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_2=East_freq_2.reset_index()
East_freq_2=East_freq_2.rename(index=str, columns={'Ind': 'IND_2' })
East_freq_3=East_freq[East_freq.Ind==3.0]
East_freq_3=East_freq_3.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_3=East_freq_3.reset_index()
East_freq_3=East_freq_3.rename(index=str, columns={'Ind': 'IND_3' })
East_freq_4=East_freq[East_freq.Ind==4.0]
East_freq_4=East_freq_4.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_4=East_freq_4.reset_index()
East_freq_4=East_freq_4.rename(index=str, columns={'Ind': 'IND_4' })
East_freq_5=East_freq[East_freq.Ind==5.0]
East_freq_5=East_freq_5.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_5=East_freq_5.reset_index()
East_freq_5=East_freq_5.rename(index=str, columns={'Ind': 'IND_5' })
East_freq_6=East_freq[East_freq.Ind==6.0]
East_freq_6=East_freq_6.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_6=East_freq_6.reset_index()
East_freq_6=East_freq_6.rename(index=str, columns={'Ind': 'IND_6' })
East_freq_7=East_freq[East_freq.Ind==7.0]
East_freq_7=East_freq_7.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_7=East_freq_7.reset_index()
East_freq_7=East_freq_7.rename(index=str, columns={'Ind': 'IND_7' })
East_freq_8=East_freq[East_freq.Ind==8.0]
East_freq_8=East_freq_8.groupby(['CHROM', 'BIN_START', 'BIN_END'])['Ind'].count()
East_freq_8=East_freq_8.reset_index()
East_freq_8=East_freq_8.rename(index=str, columns={'Ind': 'IND_8' })
windows=East_freq_1.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_2.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_3.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_4.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_5.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_6.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_7.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=East_freq_8.merge(windows, on=['CHROM', 'BIN_START', 'BIN_END'], how='outer')
windows=windows[['CHROM','BIN_START','BIN_END','IND_1','IND_2','IND_3','IND_4','IND_5','IND_6','IND_7','IND_8']]
windows=windows.fillna(0.0)


final=pandas.DataFrame(columns=['CHROM','BIN_START','Pi','Tw','Td','N_sites'])
for window_index, window_1 in windows.iterrows():
    sites_window=all_sites_2[(all_sites_2['CHROM']==window_1.CHROM)&(all_sites_2['BIN_START']==window_1.BIN_START)&(all_sites_2['BIN_END']==window_1.BIN_END)]['N_sites']
    if len(sites_window) > 0:
        sites_window=int(sites_window)
        if (sites_type=='all')&(sites_window>1000):
            sf_sum_intergenic=[window_1.IND_1,  window_1.IND_2,  window_1.IND_3,  window_1.IND_4,  window_1.IND_5,  window_1.IND_6,  window_1.IND_7,  window_1.IND_8,0,0,0,0,0,0,0, int(sites_window)]
            intergenic_stats=afs_basic_stats(sf_sum_intergenic)
            temp=pandas.DataFrame({'CHROM':[window_1.CHROM], 'BIN_START':[window_1.BIN_START],'Pi':[intergenic_stats[0]],'Tw':[intergenic_stats[1]],'Td':[intergenic_stats[2]], 'N_sites':[int(sites_window)]})
            final=final.append(temp)
        elif (sites_type=='intergenic')&(sites_window>1000):
            sf_sum_intergenic=[window_1.IND_1,  window_1.IND_2,  window_1.IND_3,  window_1.IND_4,  window_1.IND_5,  window_1.IND_6,  window_1.IND_7,  window_1.IND_8,0,0,0,0,0,0,0, int(sites_window)]
            intergenic_stats=afs_basic_stats(sf_sum_intergenic)
            temp=pandas.DataFrame({'CHROM':[window_1.CHROM], 'BIN_START':[window_1.BIN_START],'Pi':[intergenic_stats[0]],'Tw':[intergenic_stats[1]],'Td':[intergenic_stats[2]], 'N_sites':[int(sites_window)]})
            final=final.append(temp)
        elif (sites_type=='introns')&(sites_window>1000):
            sf_sum_intergenic=[window_1.IND_1,  window_1.IND_2,  window_1.IND_3,  window_1.IND_4,  window_1.IND_5,  window_1.IND_6,  window_1.IND_7,  window_1.IND_8,0,0,0,0,0,0,0, int(sites_window)]
            intergenic_stats=afs_basic_stats(sf_sum_intergenic)
            temp=pandas.DataFrame({'CHROM':[window_1.CHROM], 'BIN_START':[window_1.BIN_START],'Pi':[intergenic_stats[0]],'Tw':[intergenic_stats[1]],'Td':[intergenic_stats[2]], 'N_sites':[int(sites_window)]})
            final=final.append(temp)
        elif (sites_type in condons) & (sites_window>100) :
                sf_sum_intergenic=[window_1.IND_1,  window_1.IND_2,  window_1.IND_3,  window_1.IND_4,  window_1.IND_5,  window_1.IND_6,  window_1.IND_7,  window_1.IND_8,0,0,0,0,0,0,0, int(sites_window)]
                intergenic_stats=afs_basic_stats(sf_sum_intergenic)
                temp=pandas.DataFrame({'CHROM':[window_1.CHROM], 'BIN_START':[window_1.BIN_START],'Pi':[intergenic_stats[0]],'Tw':[intergenic_stats[1]],'Td':[intergenic_stats[2]], 'N_sites':[int(sites_window)]})
                final=final.append(temp)



base=os.path.basename(frequency_file)
base=os.path.splitext(base)[0]

final.to_csv(base+"_"+sites_type+'_'+str(part)+'_stats'+'.csv')





