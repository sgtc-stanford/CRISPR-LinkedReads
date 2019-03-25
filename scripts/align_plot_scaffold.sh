#!/bin/bash

# This script was developed and tested using bwa version: bwa-0.7.4
# User must have bwa in their path

# Function to print help
print_usage()
{
           echo "Usage: $(basename "$0") -f {fq_dir} -a  {assembly_fastq} -c {ref_contig_name}"; 
           echo "Where -f is the name of the dir containing fastq files (REQUIRED)";
	   echo "      -a is the name of the assembly fasta file (REQUIRED)";
           echo "      -c is the name of the reference contig -- can be all if want all contigs (REQUIRED)";
           echo "      -r fastq read1 (REQUIRED)";
           echo "      -s fastq read2 (REQUIRED)";
           return
}


# Parse command line options
OPTIND=1
while getopts "f:a:c:r:s:q:g:l:h" OPT
do
  case "$OPT" in
    f) fq_dir="$OPTARG";;
	a) assembly_fa="$OPTARG";;
        c) contig_ref="$OPTARG";;
        r) read1_fq="$OPTARG";;
        s) read2_fq="$OPTARG";;
        g) gen_ref="$OPTARG";;
        q) read_qual="$OPTARG";;
        l) gen_region="$OPTARG";;
    h) print_usage; exit 1;;
   \?) print_usage; exit 1;;
    :) echo "Option -$OPTARG requires an argument."; print_usage; exit 1;;
  esac
done

if [[ $# -eq 0 ]] ; then
    echo 'This script requires arguments'
    print_usage
    exit 1
fi

#gen_ref="/mnt/ix1/Resources/GenomeRef/Homo_sapiens/NCBI/GRCh38/Sequence/BWAIndex/genome.fa"

echo $fq_dir
contig_str="${contig_ref%.*}"
echo $contig_str
out_pre="$fq_dir"_"$contig_str"
out_gen="$fq_dir"_genome
echo $out_pre.bam

## Align reads to genome (if not already done)
if [ -f $out_gen.st.uq.bam ]; then
   echo "File $out_gen.st.uq.bam exists -- skipping alignment to genome."
else
	bwa mem -t 1 $gen_ref $fq_dir/$read1_fq $fq_dir/$read2_fq | samtools view -bSh - > $out_gen.bam
	samtools sort $out_gen.bam -o $out_gen.st.bam
	samtools index $out_gen.st.bam
	samtools view -b -F 256 $out_gen.st.bam -o $out_gen.st.uq.bam
	samtools index $out_gen.st.uq.bam
fi


## Align reads to the de novo contig
if [ ${assembly_fa: -3} == ".gz" ]; then
	gunzip $assembly_fa
	bgzip ${assembly_fa:: -3}
fi

#/mnt/ix1/Resources/tools/samtools-0.1.18/samtools faidx $assembly_fa
samtools faidx $assembly_fa

if [[ "$contig_str" == all ]]; then
	echo "Aligning to all contigs"
	bwa index $assembly_fa
	bwa mem -t 1 $assembly_fa $fq_dir/$read1_fq $fq_dir/$read2_fq | samtools view -bSh - > $out_pre.bam	
else
	echo "Aligning to contig $contig_str only"
	samtools faidx $assembly_fa $contig_str > $contig_str.fa
	bwa index $contig_str.fa
	## Align reads to de novo contig
	bwa mem -t 1 $contig_str.fa $fq_dir/$read1_fq $fq_dir/$read2_fq | samtools view -bSh - > $out_pre.bam
fi

samtools sort $out_pre.bam -o $out_pre.st.bam
samtools index $out_pre.st.bam
samtools view -b -F 256 $out_pre.st.bam -o $out_pre.st.uq.bam
samtools index $out_pre.st.uq.bam


## Process and merge the alignments
parse_aln_from_bam.py $out_gen.st.uq.bam $out_pre.st.uq.bam $out_pre.merge.txt

plot_scaffold.R $out_pre.merge.txt $out_pre.merge.png $gen_region $read_qual

