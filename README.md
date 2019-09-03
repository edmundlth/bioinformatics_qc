# MGHA Flagship Bioinformatics QC
----
This repository contains code to build and run a containerised environment 
to perform bioinformatics file checks for MGHA flagship data.



# Install and run instruction
----
Assume we have a directory containing genomic files (Fastq, BAM and VCF) in `data_dir` 
together with a `manifest.txt` in the same directory with at least the following columns: 
1. `sample_name`
2. `file_type` : e.g. `fastq`, `bam`, `vcf`
3. `filename`
4. `library_layout` : e.g. `Paired`

```bash
git clone https://github.com/edmundlth/bioinformatics_qc
cd bioinformatics_qc
docker build -t flagship_check:${tag} .
docker run -it --rm --volume data_dir:/mnt/ --volume <output_dir>:/home/output/ flagship_check:${tag} 
time bash run_file_check.sh
```  

