#!/usr/bin/env python

import sys
import os
import __main__ as main
import argparse
import ast

import sys
import pandas as pd
import pysam
import numpy as np

padding=50000

MIN_MAPQ = 0
PERF_CIGAR = False

bam_input = sys.argv[1]
sv_input = sys.argv[2]

bam_open = pysam.Samfile(bam_input)
sv_df = pd.read_table(sv_input, sep="\t")


def get_barcode_ids(bam_in, chrom, start, end, min_mapq, perf_cigar):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq: #and (not(perf_cigar) or (not(r.cigar is None) and len(r.cigar) == 1)):
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return list(set(list(bcs)))


df_list_bcs = []
df_list_nobcs = []

for index, row in sv_df.iterrows():
	print row['name']

	# initialize lists to store output
	bam_in_regions = []
	bam_out_regions = []

	# get info from the data frame
	sv_name,subname,chrom1,start1,stop1,chrom2,start2,stop2,sv_type,bc_filt=str(row['name']),str(row['subname']),str(row['chrom1']),int(row['start1']),int(row['stop1']),str(row['chrom2']),int(row['start2']),int(row['stop2']),str(row['type']),str(row['bc_filt'])

	# depending on the SV type, define the "in" and "out" regions
	if sv_type=="DEL":
		bam_in_regions.append(";".join([chrom1,str(start1),str(stop1)]))	
		bam_in_regions.append(";".join([chrom2,str(start2),str(stop2)]))
		bam_out_regions.append(";".join([chrom1,str(stop1+200),str(start2-200)]))

	elif sv_type=="INV_LEFT":
		window_1 = stop1-start1
		window_1_adj = min(window_1,10000)
		window_2 = stop2-start2
		window_2_adj = min(window_2,10000)		
		
		bam_in_regions.append(";".join([chrom1,str(stop1-window_1_adj),str(stop1)]))	
		bam_in_regions.append(";".join([chrom2,str(stop2-window_2_adj),str(stop2)]))
		bam_out_regions.append(";".join([chrom2,str(stop2+200),str(stop2+50000)]))
		
	elif sv_type=="INV_RIGHT":
		window_1 = stop1-start1
		window_1_adj = min(window_1,10000)
		window_2 = stop2-start2
		window_2_adj = min(window_2,10000)			
	
		bam_in_regions.append(";".join([chrom1,str(start1),str(start1+window_1_adj)]))	
		bam_in_regions.append(";".join([chrom2,str(start2),str(start2+window_2_adj)]))
		bam_out_regions.append(";".join([chrom1,str(start1-50000),str(start1-200)]))
	
	else:
		print "problem: dont't handle sv_type: " + str(sv_type)

	# get barcodes inside each "in" region
	bam_in_bcs=[]
	for ri in bam_in_regions:
		chr = str(ri.split(';')[0])
		start = int(ri.split(';')[1])
		end = int(ri.split(';')[2])
		region_bc_list = get_barcode_ids(bam_open, chr, start, end, MIN_MAPQ, PERF_CIGAR)
		bam_in_bcs.append(region_bc_list)

	# get barcodes shared by all "in" regions
	bam_in_bcs_comm = set(bam_in_bcs[0])
	for s in bam_in_bcs[1:]:
		bam_in_bcs_comm.intersection_update(s)

	# get barcodes inside each "out" region
	bam_out_bcs=[]
	for ro in bam_out_regions:
		chr = str(ro.split(';')[0])
		start = int(ro.split(';')[1])
		end = int(ro.split(';')[2])
		region_bc_list = get_barcode_ids(bam_open, chr, start, end, MIN_MAPQ, PERF_CIGAR)
		bam_out_bcs.append(region_bc_list)
	
	# flatten list of bcs in "out" regions
	bam_out_bcs = [val for sublist in bam_out_bcs for val in sublist]
	#bam_out_bcs = [b.split("-")[0] for b in bam_out_bcs]

	# find those barcodes in the "in" regions and NOT in the "out" regions
	all_bcs = list(set(bam_in_bcs_comm) - set(bam_out_bcs))
	#all_bcs = [b.split("-")[0] for b in all_bcs]
	
	# add data to lists for output
	if bc_filt == "in_out":
		df_list_bcs.append([chrom1,start1-padding,stop2+padding,subname,tuple(all_bcs)])
	elif bc_filt == "out":
		df_list_bcs.append([chrom1,start1-padding,stop2+padding,subname,tuple(bam_out_bcs)])
	elif bc_filt == "in":
		df_list_bcs.append([chrom1,start1-padding,stop2+padding,subname,tuple(bam_in_bcs_comm)])
	else:
		print "bc_filt columns must contain one of: in_out,out,in"
		sys.exit()
	
	#df_list_bcs.append([chrom1,start1,stop1,chrom2,start2,stop2,sv_type,bam_in_bcs_comm,bam_out_bcs,all_bcs,len(bam_in_bcs_comm),len(bam_out_bcs),len(all_bcs), row['name'],row['sv_name']])
	#df_list_nobcs.append([chrom1,start1,stop1,chrom2,start2,stop2,sv_type,len(bam_in_bcs_comm),len(bam_out_bcs),len(all_bcs), row['name'],row['sv_name']])

# turn output lists into data frames
df_out_bcs = pd.DataFrame(df_list_bcs, columns=['chr','start','end','name','bcs'])
#df_out_nobcs = pd.DataFrame(df_list_nobcs, columns=['chrom1','start1','stop1','chrom2','start2','stop2','sv_type','num_bam_in_bcs_comm','num_bam_out_bcs','num_all_bcs','name','sv_name'])

# write output dfs to files
df_out_bcs.to_csv("shared_bcs.txt", sep="\t", index=False)
#df_out_nobcs.to_csv("shared_nobcs.txt", sep="\t", index=False)

