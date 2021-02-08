import sys
import os
import string
import pandas
import numpy as np
import matplotlib.pyplot as plt
import statsmodels.api as sm
import seaborn as sns
import itertools
import math
import time
import sys
import personal_popgen
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqUtils import GC
import scipy


#get genes that are in one orf, start to stop codon, 

fasta_DPlex=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ALL_induveduals.CDS.fasta')
fasta_DErip=personal_popgen.fasta_dict('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/Outgroups/danaus_eresimus.CDS.fasta')

good_genes=[]
realign=[]
###get gene wize allignments
for i in fasta_DPlex.keys():
	if i in (fasta_DErip.keys()) and (len(fasta_DErip[i]) > 0.0):
		Ns=float(str(fasta_DErip[i]).count('N'))/len(fasta_DErip[i])
		if Ns < 0.35:
			if len(fasta_DPlex[i]) != len(fasta_DErip[i]):
				realign.append(i[0:-3])
			elif len(fasta_DPlex[i]) == len(fasta_DErip[i]):
				good_genes.append(i[0:-3])
			#F = open('/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/dnds/Danaus_eresimus/gene_allignments/'+i+'.CDS.fasta',"w") 
			#F.write(">"+str('Danaus_plexippus')+"\n")
			#F.write(str(fasta_DPlex[i])+"\n")
			#F.write(">"+str('Danaus_eresimus')+"\n")
			#F.write(str(fasta_DErip[i])+"\n")
			#F.close()

####load gene details

head_stuff=['scaffold','source','feature','start','end','.','orentation','..','infor']
annotation=pandas.read_table("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/genome/Danaus_plexippus_v3_-_genes.gff",skiprows=0, names=head_stuff)
annotation=annotation[['scaffold','source','feature','start','end','infor','orentation']]
annotation_gene=annotation[annotation.feature=='gene']
annotation_gene['ID']=annotation_gene.infor.str.split(';').str.get(0).str.split('=').str.get(1)

annotation_gene=annotation_gene[['scaffold', 'start', 'end', 'orentation','ID']]
annotation_gene=annotation_gene.rename(columns={'scaffold': 'CHROM', 'start': 'BIN_START', 'end': 'BIN_END'})
annotation_gene['Sites']=0
annotation_gene['Sites']=annotation_gene['BIN_END']-annotation_gene['BIN_START']
annotation_gene=annotation_gene[['CHROM','BIN_START','BIN_END','Sites','ID']]

####load gene pn/ps details per gene

fianl_df_5=pandas.read_csv( '/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/temp/ALL_induveduals_pn_ps.pnps',sep=" ",index_col=False)[["geneID","Pin","Pis","length_orf","Pi_syn","Pi_nonsyn","syn_sites","nonsyn_sites"]]
fianl_df_5=fianl_df_5.rename(columns={'geneID': 'Gene_ID'})
fianl_df_5['Gene_ID']=fianl_df_5['Gene_ID'].str[0:-3]
fianl_df_5=fianl_df_5.rename(columns={ 'Pin': 'Pin_Monarch', 'Pis': 'Pis_Monarch', "Pi_nonsyn": "Pi_nonsyn_Monarch", "Pi_syn": "Pi_syn_Monarch","nonsyn_sites":"nonsyn_sites_Monarch", "syn_sites": "syn_sites_Monarch"})
fianl_df_5['Pin_Pis_Monarch']=fianl_df_5["Pin_Monarch"]/fianl_df_5["Pis_Monarch"]
fianl_df_5['Pin_Pis_Monarch']=fianl_df_5['Pin_Pis_Monarch'].replace(np.inf, np.nan)

####load gene dn/ds details per gene

dnds_table=pandas.read_csv("/proj/uppstore2017185/b2014034_nobackup/Venkat/Monarch_stuff/dnds/Danaus_eresimus/gene_allignments/East/dn_ds_1/final_dnds_result",  sep='\t')

dnds_table['W_Danaus_plexippus']=dnds_table['dn']/dnds_table['ds']
dnds_table['s_subs']=dnds_table['ds']*dnds_table['N_sites']
dnds_table['n_subs']=dnds_table['dn']*dnds_table['S_sites']

dnds_table['geneID']=dnds_table.gene_name.str.split('.').str.get(-3).str[1:-3]
dnds_table=dnds_table[['geneID','dn','ds','N_sites','S_sites','W_Danaus_plexippus','s_subs','n_subs']]
dnds_table=dnds_table.rename(columns={'geneID': 'Gene_ID'})
dnds_table['Gene_ID']=dnds_table['Gene_ID']
dnds_table['W_Danaus_plexippus']=dnds_table['W_Danaus_plexippus'].replace(np.inf, np.nan)

dnds_table=dnds_table.merge(fianl_df_5, on = ['Gene_ID'])
dnds_table['DOS']=(dnds_table['dn']/(dnds_table['dn']+dnds_table['ds'])) - (dnds_table['Pin_Monarch']/(dnds_table['Pin_Monarch']+dnds_table['Pis_Monarch']))

dnds_table['Pin_Pis_Monarch']=dnds_table['Pin_Pis_Monarch'].replace(np.inf, np.nan)


dnds_table['alpha']=1-(dnds_table['Pin_Pis_Monarch']/dnds_table['W_Danaus_plexippus'])
dnds_table['W_alpha']=dnds_table['alpha']*dnds_table['W_Danaus_plexippus']

dnds_table['alpha']=dnds_table['alpha'].replace(np.inf, np.nan)
dnds_table['alpha']=dnds_table['alpha'].replace(-np.inf, np.nan)
dnds_table['W_alpha']=dnds_table['W_alpha'].replace(np.inf, np.nan)
dnds_table['W_alpha']=dnds_table['W_alpha'].replace(-np.inf, np.nan)


