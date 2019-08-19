#!/usr/bin/env bash

# Get input arguments
sam_file=$1
vcf_file=$2
#output_file=$3

# Grab samtools header
sam_file_contigs (){
samtools view -H "${1}" | {
    # Collect the SQ tag
    grep '^@SQ'
} | {
    # Grab only columns two and three
    cut -f 2,3
} | {
    # Remove the ln: and sq: attributes
    sed 's/\wN://g'
} | {
    # Sort chromosomes alphabetically and write to output sam file
    sort -k1
}
}

vcf_file_contigs (){
# Generate contigs and lengths from vcf
cat "${1}" |  {
    # Get contigs header
    grep '^##contig'
} | {
    # Strip ID prefix
    sed 's/##contig=<ID=//'
} | {
    # Strip length prefix
    sed 's/length=//'
} | {
    # Get columns 1 and 2
    cut -d',' -f1,2
} | {
    # Convert commas to tabs
    tr ',' '\t'
} | {
    # Sort chromosomes alphabetically and write to output sam file
    sort -k1
}
}



sam_contigs=$(sam_file_contigs ${sam_file})
vcf_contigs=$(vcf_file_contigs ${vcf_file})
if [ "$sam_contigs" == "$vcf_contigs" ]
then
    exit 0
else
    exit 1
fi