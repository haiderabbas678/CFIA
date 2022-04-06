# Genome_Assembly
# Introduction
This script is able to assemble genomes(de novo) using illumina paired end reads as well as produce a QC report to assess the quality of the assembled /
genome

# Dependencies
conda install -c bioconda fastp=0.23.2/
conda install -c bioconda skesa=2.4.0/
conda install -c bioconda bwa=0.7.17/
conda install -c bioconda samtools= 1.15/
conda install -c bioconda qualimap=2.2.2d

# How to execuet the script in terminal
Go into your desired directory and execute the below command 
python genome_assembly.py [fold]/
[fold] The folder that conatins the query files of interest

