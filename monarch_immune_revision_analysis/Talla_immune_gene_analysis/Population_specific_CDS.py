import sys
import os
import string
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
import Bio.Data.CodonTable
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq
import argparse
import pandas


def join_raw_data_base(lists):
    for i in range(1,len(lists)):
        j=i-1
        if i == 1:
            new_2=pandas.merge(lists[j], lists[i], on=['CHROM','BIN_START','BIN_END'], how='inner', suffixes=('_1', '_2'))
        else:
            new_2=pandas.merge(nextone, lists[i], on=['CHROM','BIN_START','BIN_END'], how='inner')

        nextone=new_2
    return nextone

def join_CDS_seq(lists):
    seq=''
    if lists[0] == '-':
        list_to_use = lists[1:][::-1]
        for i in list_to_use:
            my_seq = Seq(str(i))
            seq += str(my_seq.reverse_complement()) 
    elif lists[0] == '+':
        list_to_use = lists[1:]
        for i in list_to_use:
            seq += str(i)
    return seq

def get_intron(CDS):
    df1=CDS
    df1.index = range(1, 2*len(df1)+1, 2)
    df2 = df1.reindex(index=range(2*len(df1)))
    df2 = df2.iloc[1:]
    df2.scaffold = df2.scaffold.ffill()
    df2.source = df2.source.ffill()
    df2.feature = df2.feature.fillna('intron')
    df2.infor = df2.infor.ffill()
    int_starts=list(df2.end.dropna())[:-1]
    int_end=list(df2.start.dropna())[1:]
    df2.loc[df2.start.isnull(), 'start'] = int_starts
    df2.loc[df2.end.isnull(), 'end'] = int_end
    return df2


def get_args():
    parser = argparse.ArgumentParser(description='''Outputs CDS for the population from allele freq file

        -f freq file
        -i number of chromosomes (induvuduals in the population *2)
        -o Output
        ''')
    parser.add_argument('-f', '--freq', type=str, required=True)
    parser.add_argument('-o', '--out', type=str, required=True)
    parser.add_argument('-i', '--chr', type=int, required=True)
    return parser.parse_args()


def output_good_sites_fequency(freq_file1,chr_filter):
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
    #temp2=temp2[temp2.minor_freq>0]
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
    #temp2=temp2[temp2.minor_freq>0]
    temp2=temp2[['CHROM', 'POS', 'major_allele', 'minor_allele','major_freq', 'minor_freq']]
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
    return new_output



args = get_args()


sns.set(font_scale=1.5)
sns.set_style("whitegrid", {'axes.grid' : False})



CDS_fasta=personal_popgen.fasta_dict('/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_cds.fa')
goods=output_good_sites_fequency(args.freq,args.chr)
#goods=output_good_sites_fequency("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ALL_induveduals.frq",86)
goods=goods[(goods.major_allele!='*') & (goods.minor_allele!='*')]
goods['POS']=goods['POS'].astype(int)
fixed=goods[(goods['minor_freq']==0.0)]
fixed=fixed.append(goods[(goods['minor_freq']>=round(float(2.0/(args.chr)), 8 ))], ignore_index = True)
goods=fixed
goods['CODE']=''
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'AG') | (goods['minor_allele']+goods['major_allele'] == 'GA')) ]='R'
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'CT') | (goods['minor_allele']+goods['major_allele'] == 'TC')) ]='Y'
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'GC') | (goods['minor_allele']+goods['major_allele'] == 'CG')) ]='S'
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'AT') | (goods['minor_allele']+goods['major_allele'] == 'TA')) ]='W'
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'GT') | (goods['minor_allele']+goods['major_allele'] == 'TG')) ]='K'
goods['CODE'][ (goods['minor_freq'] >= round(float(2.0/(args.chr)), 8 )) & ((goods['minor_allele']+goods['major_allele'] == 'AC') | (goods['minor_allele']+goods['major_allele'] == 'CA')) ]='M'



goods['CODE'][ (goods['major_freq'] == 1.00  )]=goods[ (goods['major_freq'] == 1.00  )]['major_allele']

fasta=personal_popgen.fasta_dict('/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_scaffolds.fa')

head_stuff=['scaffold','source','feature','start','end','extra','-','indexing','infor']
#CDSs=pandas.read_table("/home/venkat/bin/snpgenie/test/",skiprows=1, names=head_stuff)

annotation=pandas.read_table("/work/smalab/Venkat/genome/Danaus_plexippus_v3_-_genes.gff",skiprows=1, names=head_stuff)
annotation=annotation[['scaffold','source','feature','start','end','-','indexing','infor']]

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='+'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='+'), 'start']-1

annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='1') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='2') & (annotation['-']=='-'), 'start']-1

annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='-'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='-'), 'start']-1
annotation.loc[(annotation['indexing']=='0')  & (annotation['-']=='+'), 'start']=annotation.loc[(annotation['indexing']=='0') & (annotation['-']=='+'), 'start']-1

#annotation_gene=annotation[(annotation.feature=='gene') & (annotation.infor.str[31:37]==';Name=')]
annotation_gene=annotation[(annotation.feature=='gene')]
#annotation_gene=pandas.read_csv(str(sys.argv[1]))


CDS_cordinates=pandas.DataFrame(columns=['start','end','gene_name'])

#F = open('ALL_induveduals'+'.CDS.fasta',"w") 
F = open(args.out+'.CDS.fasta',"w") 

for scaffold in set(list(annotation_gene['scaffold'])):
    sequence=MutableSeq(str(fasta[scaffold]))
    scaffold_gene=annotation_gene[(annotation_gene.scaffold == scaffold)]
    for gene_index, gene in scaffold_gene.iterrows():
        marker_name=gene.infor.split(';')[-1].split('=')[-1]
        gene_all_temp=annotation[(annotation.start >= int(gene['start'])) & (annotation.end <= int(gene['end'])) & (gene['scaffold'] == annotation.scaffold)]
        gene_all=annotation[(annotation.infor.str.split(';').str.get(-1).str.split('=').str.get(-1).str.match(marker_name))]
        CDS_name=list(gene_all_temp[gene_all_temp.feature=='mRNA'].infor.str.split(';').str.get(0).str.split('=').str.get(1))[0]
        CDS=gene_all[(gene_all.feature=='CDS')]
        exon_stuff=gene_all[(gene_all.feature=='exon')]
        CDS=CDS[(CDS.infor.str[3:len(CDS_name)+3]==CDS_name )]
        exon_stuff=get_intron(exon_stuff)
        exon_stuff.loc[exon_stuff.feature=='intron','start']=exon_stuff.loc[exon_stuff.feature=='intron','start']+1
        exon_stuff.loc[exon_stuff.feature=='intron','end']=exon_stuff.loc[exon_stuff.feature=='intron','end']-1
        #CDS=CDS.append(exon_stuff[exon_stuff.feature=='intron'], ignore_index=True)
        gene_sequnce_cds=CDS_fasta[CDS_name]
        gene_sequnce_extract_list=[gene['-']]
        CDS['diff']=CDS['end']-CDS['start']
        CDS['CDS_start']=0
        CDS['CDS_end']=0
        if len(CDS) > 0:
            if CDS['-'].unique()[0]=='-':
                CDS=CDS.reindex(index=CDS.index[::-1])
                CDS=CDS.reset_index()
                CDS_start=[0]
                CDS_end=[]
                for i in list(CDS['diff']):
                    CDS_start.append(CDS_start[-1]+i)
                CDS_end=CDS_end+CDS_start[1:]
                CDS_start=CDS_start[:-1]
                CDS['CDS_start']=pandas.Series(CDS_start)
                CDS['CDS_end']=pandas.Series(CDS_end)
            elif CDS['-'].unique()[0]=='+':
                CDS=CDS.reset_index()
                CDS_start=[0]
                CDS_end=[]
                for i in list(CDS['diff']):
                    CDS_start.append(CDS_start[-1]+i)
                CDS_end=CDS_end+CDS_start[1:]
                CDS_start=CDS_start[:-1]
                CDS['CDS_start']=pandas.Series(CDS_start)
                CDS['CDS_end']=pandas.Series(CDS_end)
            #print CDS
        sequence_final=''
        for n_i,n in CDS.iterrows():
            if n['-'] == '-':
                change=goods[(goods['CHROM']==n.scaffold) & (goods['POS']>=n.start) & (goods['POS']<=n.end)]
                if len(change) > 0:
                    for m_i, m in change.iterrows():
                        sequence[int(m.POS)-1]=m['CODE']
                sequence_final+=str(Seq(str(sequence[n.start:n.end])).reverse_complement())
            elif n['-'] == '+':
                change=goods[(goods['CHROM']==n.scaffold) & (goods['POS']>=n.start) & (goods['POS']<=n.end)]
                if len(change) > 0:
                    for m_i, m in change.iterrows():
                        sequence[int(m.POS)-1]=m['CODE']
                sequence_final+=str(sequence[n.start:n.end])
        F.write(">"+str(CDS_name)+"\n")
        F.write(str(sequence_final)+"\n")
        sys.stdout.flush()
        #print ">"+CDS_name
        #print sequence_final

