#!/usr/bin/env python

import sys, pysam, pandas as pd

infile_ref = sys.argv[1]
infile_asm = sys.argv[2]
outfile = sys.argv[3]

bamfile_ref = pysam.AlignmentFile(infile_ref, "rb")
bamfile_asm = pysam.AlignmentFile(infile_asm, "rb")

read_info_ref = []

for read in bamfile_ref.fetch():
	cig_str = read.cigarstring
	cig_tup = read.cigar
	ref_id=bamfile_ref.get_reference_name(read.reference_id)
	read_info_ref.append(eval('read.query_name,ref_id,read.reference_start+1,len(read.query_sequence),read.mapping_quality,cig_str,cig_tup'))
	
bamfile_ref.close()

df_ref = pd.DataFrame(read_info_ref, columns=['read','chr','pos','len','qual','cigstr','cigtup'])


read_info_asm = []

for read in bamfile_asm.fetch():
        cig_str = read.cigarstring
        cig_tup = read.cigar
        ref_id=bamfile_asm.get_reference_name(read.reference_id)
        read_info_asm.append(eval('read.query_name,ref_id,read.reference_start+1,len(read.query_sequence),read.mapping_quality,cig_str,cig_tup'))
        
bamfile_asm.close()

df_asm = pd.DataFrame(read_info_asm, columns=['read','chr','pos','len','qual','cigstr','cigtup'])


## merge below

df = pd.merge(df_ref, df_asm, on="read", how="inner")
df.to_csv(outfile, sep="\t", index=False)

