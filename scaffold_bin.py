#!/usr/bin/env python

import sys
import pandas as pd
import re
import pysam
import numpy as np
from collections import OrderedDict
import faulthandler
import csv
import mappy as mp
from Bio import SeqIO
import gzip

infile = sys.argv[1]
genome = sys.argv[2]
phased_basic_file = sys.argv[3]
filt_len = sys.argv[4]
filt_n = sys.argv[5]


# alignment file
#infile="MHC_4_asmbly.filt.fasta"
#genome="/mnt/ix1/Resources/GenomeRef/Homo_sapiens/NCBI/GRCh38/Sequence/MinimapIndex/genome.mmi"
outfile_list = infile.split(".")[:-1]
outfile = "_".join(outfile_list) + ".bin.txt"

asm_preset = "asm5"
outfile_full_name = "_".join(outfile_list) + ".tmp.txt"

# phased basic file
#phased_basic_file = "hg38.hybrid.phased_basic.gz"
tbx_file = pysam.TabixFile(phased_basic_file)

a = mp.Aligner(str(genome), preset=asm_preset)

if not a: raise Exception("ERROR: failed to load/build index")

outfile_full = open(outfile_full_name, 'w')
writer = csv.writer(outfile_full, delimiter='\t')
writer.writerows([['contig','chr','genome_start','genome_end','mapq','is_primary','block_id','genome_pos','hap1_compare','hap2_compare','hap1_ppn','hap2_ppn']])

for h in tbx_file.header:
	tbx_header=h.split('\t')
	tbx_header=[x.replace('#','') for x in tbx_header]

def check_match(t):
	genome_s = str(t['genome_seq'])
	contig_s = str(t['contig_seq'])
	
	if (genome_s=="-" or contig_s=="-"):
		return "-"
	elif genome_s==contig_s:
		return "+"
	elif genome_s!=contig_s:
		return "*"
	else:
		return "problem"

def get_tbx_info(r,h):
	
	tbx_row_list = r.split('\t')
	tbx_row_dict = dict(zip(h, tbx_row_list))
	tbx_chrom = str(tbx_row_dict.get('chrom'))
	tbx_pos = int(tbx_row_dict.get('pos'))
	tbx_id = str(tbx_row_dict.get('block_id'))
	tbx_base1 = str(tbx_row_dict.get('base_1'))
	tbx_base2 = str(tbx_row_dict.get('base_2'))
	tbx_var = str(tbx_row_dict.get('var_type'))
	tbx_phase = str(tbx_row_dict.get('phase_status'))
	tbx_hom = str(tbx_row_dict.get('hom_status'))
	tbx_filt = str(tbx_row_dict.get('filter'))
	tbx_ref = str(tbx_row_dict.get('ref'))

	#if (tbx_var=='snv' and tbx_filt=='[]' and ((tbx_phase=='phased' and tbx_hom=='het') or (tbx_hom=='hom'))):
	if (tbx_var=='snv' and tbx_filt=='[]' and (tbx_phase=='phased' and tbx_hom=='het')):
		return [tbx_chrom, tbx_pos, tbx_id, tbx_ref, tbx_base1, tbx_base2]

#def parse_cs(row,tbx):

## filter the input fasta file

keep_seqs = []

if infile.endswith(".gz"):
	outfile_list = infile.split(".")[:-2]
	outfile_fa = "_".join(outfile_list) + ".filt.fasta"
	print outfile
	with gzip.open(infile, "rt") as handle:
		for rec in SeqIO.parse(handle, 'fasta'):
			name = rec.id
			seq = rec.seq
			seqLen = len(rec)
			if int(seqLen)>int(filt_len):
				nCount = seq.upper().count('N')
				perc_n = (float(nCount)/seqLen)*100
				if perc_n < float(filt_n):
					print name, seqLen, nCount, perc_n
					keep_seqs.append(rec)

else:
	outfile_list = infile.split(".")[:-1]
	outfile_fa = "_".join(outfile_list) + ".filt.fasta"
	FastaFile = open(infile, 'rU')
	for rec in SeqIO.parse(FastaFile, 'fasta'):
		name = rec.id
		seq = rec.seq
		seqLen = len(rec)
		if int(seqLen)>int(filt_len):
			nCount = seq.upper().count('N')
			perc_n = (float(nCount)/seqLen)*100
			if perc_n < float(filt_n):
				print name, seqLen, nCount, perc_n
				keep_seqs.append(rec)
	FastaFile.close()


SeqIO.write(keep_seqs, outfile_fa, "fasta")

## end of filtering ##

## now read in and process the 'filtered' fasta files

for name, seq, qual in mp.fastx_read(outfile_fa):
	print name
	seqs = re.split('N+', seq) 

	z=0
	for s in seqs:
		
		z = z+1
		
		new_name = str(name) + "_" + str(z)
		
		for hit in a.map(s, cs=True):

			chr = str(hit.ctg)

			if chr in ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']:

				read = new_name
				cs =  str(hit.cs).upper()
				r_st = int(hit.r_st)
				r_en = int(hit.r_en)
				q_st = int(hit.q_st)
				seq = str(s)

				print read

				i = q_st #initialize contig start pos
				m = r_st+1 #initialize genome start pos

				## Create a data frame that summarizes cs for each alignment -- goal is to recreate ref and contig seq in alignment region

				cs_split = re.findall('[=|:|\*|\+|\-|~][a-zA-Z0-9]*',cs)
	
				chr_list_1 = []
				c_pos_list_1 = []
				c_seq_list_1 = []
				g_pos_list_1 = []
				g_seq_list_1 = []
	
				for c in cs_split:
		
					c_op=c[0] #operation
					c_det=c[1:] #details

					if c_op==":": #identical sequence
			
						j=i+int(c_det) 
						g_seq = list(seq[i:j])
						c_seq = list(seq[i:j])
						n=m+int(c_det)
			
						align_len = n-m
						chr_list = [str(chr)] * (align_len)
						contig_pos = range(i,j)
						genome_pos = range(m,n)
			
						chr_list_1=chr_list_1+chr_list
						c_pos_list_1=c_pos_list_1+contig_pos
						c_seq_list_1=c_seq_list_1+c_seq
						g_pos_list_1=g_pos_list_1+genome_pos
						g_seq_list_1=g_seq_list_1+g_seq

						i=j
						m=n
			
					if c_op=="*": #substitution
		
						g_seq=list(c_det[0])
						c_seq=list(c_det[1])
						n=m+1
						j=i+1

						align_len = n-m
						chr_list = [str(chr)] * (align_len)
						contig_pos = range(i,j)
						genome_pos = range(m,n)
					
						chr_list_1=chr_list_1+chr_list
						c_pos_list_1=c_pos_list_1+contig_pos
						c_seq_list_1=c_seq_list_1+c_seq
						g_pos_list_1=g_pos_list_1+genome_pos
						g_seq_list_1=g_seq_list_1+g_seq
						i=j
						m=n
			
					if c_op=="-": #insertion to ref

						g_ins_len=len(c_det)
						g_seq=list(str(c_det))
						c_seq=list('-' * g_ins_len)
						n=m+g_ins_len

						align_len = n-m
						chr_list = [str(chr)] * (align_len)
						contig_pos = list('-' * align_len)
						genome_pos = range(m,n)
					
						chr_list_1=chr_list_1+chr_list
						c_pos_list_1=c_pos_list_1+contig_pos
						c_seq_list_1=c_seq_list_1+c_seq
						g_pos_list_1=g_pos_list_1+genome_pos
						g_seq_list_1=g_seq_list_1+g_seq
						m=n
			
					if c_op=="+": #deletion from ref
		
						g_del_len=len(c_det)
						c_seq=list(str(c_det))
						g_seq=list('-' * g_del_len)
						j=i+g_del_len

						align_len = j-i
						chr_list = list('-' * align_len)
						contig_pos = range(i,j)
						genome_pos = list('-' * align_len)

						chr_list_1=chr_list_1+chr_list
						c_pos_list_1=c_pos_list_1+contig_pos
						c_seq_list_1=c_seq_list_1+c_seq
						g_pos_list_1=g_pos_list_1+genome_pos
						g_seq_list_1=g_seq_list_1+g_seq
						i=j
		
					if c_op=="~":
			
						print "This function does not currently support ~"
						sys.exit()

				df_final_dict = {'genome_chr': chr_list_1, 'genome_pos':g_pos_list_1,'genome_seq':g_seq_list_1,'contig_pos':c_pos_list_1,'contig_seq':c_seq_list_1}
				df_final = pd.DataFrame(data=df_final_dict)
				df_final = df_final[['genome_chr','genome_pos','genome_seq','contig_pos','contig_seq']]
				df_final['check_match'] = df_final.apply(lambda s: check_match(s), axis=1)
				df_final = df_final.loc[df_final['check_match']!="-"]
				df_final = df_final.loc[df_final['contig_seq']!="N"]
	
				# Here, bring in tabix
				tbx_list = []
				for tbx_row in tbx_file.fetch(str(chr), int(r_st), int(r_en)):
					tbx_info = get_tbx_info(tbx_row, tbx_header)
					if str(tbx_info)!="None":
						tbx_list.append(tbx_info)
	
				df_tbx = pd.DataFrame(tbx_list, columns=['genome_chr','genome_pos','block_id','ref','base_1','base_2'])
	
				# Merge the data!!
	
				df_tbx[['genome_chr']] = df_tbx[['genome_chr']].astype(str)
				df_tbx[['genome_pos']] = df_tbx[['genome_pos']].astype(int)
			
				df_final[['genome_chr']] = df_final[['genome_chr']].astype(str)
				df_final[['genome_pos']] = df_final[['genome_pos']].astype(int)
				df_merge = pd.merge(df_tbx, df_final, on=['genome_chr','genome_pos'],how='left')
	
				if not df_merge.empty:
	
					df_merge['hap1_compare'] = df_merge.apply(lambda v: 1 if v['contig_seq']==v['base_1'] else 0, axis=1)
					df_merge['hap2_compare'] = df_merge.apply(lambda v: 1 if v['contig_seq']==v['base_2'] else 0, axis=1)
					grouped = df_merge.groupby(['genome_chr','block_id']).agg(OrderedDict([('genome_pos', ['min','max','count']) , ('hap1_compare' , np.sum),('hap2_compare' , np.sum)])).reset_index()
					grouped.columns = ["_".join(x) for x in grouped.columns.ravel()]
					grouped.columns = ['chr','block_id','genome_start','genome_end','genome_pos','hap1_compare','hap2_compare']
					grouped['contig'] = read
					grouped['hap1_ppn'] = grouped.apply(lambda y: float(y['hap1_compare'])/float(y['genome_pos']) if float(y['genome_pos'])>0 else 'na', axis=1)
					grouped['hap2_ppn'] = grouped.apply(lambda y: float(y['hap2_compare'])/float(y['genome_pos']) if float(y['genome_pos'])>0 else 'na', axis=1)
					grouped['mapq'] = hit.mapq
					grouped['is_primary'] = hit.is_primary
					grouped = grouped[['contig','chr','genome_start','genome_end','mapq','is_primary','block_id','genome_pos','hap1_compare','hap2_compare','hap1_ppn','hap2_ppn']]
					grouped_list = grouped.values.tolist()
					writer.writerows(grouped_list)

				else:
					grouped_list = 	[[read,chr,r_st,r_en,hit.mapq,hit.is_primary,'na','na','na','na','na','na']]
					writer.writerows(grouped_list)

outfile_full.close()

# read back in the file  generated above and summarize the results by contig

df = pd.read_table(outfile_full_name, sep="\t")

df['scaffold'] = df['contig'].apply(lambda x: x.split("_")[0])
df[['is_primary']] = df[['is_primary']].astype(str)
df_summ = df.loc[(df['is_primary']=="True") & (df['genome_pos']!="na")]

df_summ[['genome_pos','hap1_compare','hap2_compare']] = df_summ[['genome_pos','hap1_compare','hap2_compare']].astype(int)
df_summ = df_summ[['scaffold','genome_pos','hap1_compare','hap2_compare']]
df_grp = df_summ.groupby('scaffold').agg({'genome_pos': 'sum', 'hap1_compare': 'sum', 'hap2_compare': 'sum'}).reset_index()
df_grp['num_shared_pos'] = df_grp['hap1_compare'] + df_grp['hap2_compare'] 
df_grp['hap1_perc'] = df_grp.apply(lambda row: row['hap1_compare']/float(row['num_shared_pos'])*100, axis=1)
df_grp['hap2_perc'] = df_grp.apply(lambda row: row['hap2_compare']/float(row['num_shared_pos'])*100, axis=1)
df_grp['hap_col'] = df_grp[['hap1_compare','hap2_compare']].idxmax(axis=1)
df_grp['hap_assignment'] = df_grp['hap_col'].apply(lambda x: x.split("_")[0])
df_grp = df_grp[['scaffold','genome_pos','num_shared_pos','hap1_compare','hap2_compare','hap1_perc','hap2_perc','hap_assignment']]

df_grp.to_csv(outfile, sep="\t", index=False)


