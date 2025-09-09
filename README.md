# manzanaa
Mutation Accumulation, Network, and Zoning Analysis for Nucleotides and Amino Acids

Welcome to MANZANAA:
Mutation Accumulation, Network, and Zoning Analysis for Nucleotides and Amino Acids

The following is a guide for utilizing the scripts within this project. It is highly recommended to set up all packages required to operate the scripts within a virtual environment or package manager (e.g. anaconda, miniconda, and the like).

There are two main directories within the project root, titled "manzanaa_analysis" and "manzanaa_data", respectively. The "analysis" folder contains the scripts to be executed. When you're ready to execute, the "src" folder within "analysis" is the directory to which should be navigated in your terminal.
    Note: for paired end Illumina data, use the "01_" script, which includes a module to merge "R1" and "R2" reads from the same sample. for nanopore-like data (i.e. long read sequencing), use the "01b_" script, which withholds the merging function but performs all other relevant operations.

The "data" directory is split into two sub-directories: "manzanaa_outputs" and "ngs_raw". Please format your "ngs_raw" by experiment iteration sub-folders (historically in the format "fl_XXX", with first and last initial as the two letter prefix and a 3 digit experiment iterator -- if you reach more than 999 jobs, I commend you).

Within each experiment subfolder, create the following directories:
    "demultiplex": place your fastq files (gzipped) here and a fl_xxx_demux_metadata.txt (see example metadata and separate metadata readme)
    "references": place all relevant reference .fasta files (will be called on for each entry in the metadata file)

Within "manzanaa_outputs" the scripts will generate an experiment iteration sub-folder matching the directory name in "ngs_raw". The following outputs will be populated within subfolders:
    1. fastqc_output: fastqc reports for each sample fastq (both zipped raw data and .html file)
    2. fastp_output: fastp .html report, .json output. and a filtered fastq removing reads with low quality bases
    3. merged_reads: output of NGmerge combining R1 and R2 paired end reads. Will also output reads that do not pass NGmerge as a separate .fastq
    4. seqkit_stats: filtering statistics on number of reads per barcode (if in use)
    5. alignments: holds .bam files generated from minimap2 which will be converted to dataframes in downstream scripts.
    6. unmapped: within the alignments folder. fastqc analysis on the unmapped reads that did not pass through minimap2.

for script "01b_", all above outputs will be generated with the exception of "merged_reads", as there are no paired reads to merge.

Both scripts have several tunable parameters within each module (e.g. NGmerge, minimap2). It is encouraged to adjust the desired parameters as guided by each module's documentation.

script "02_" is a python script made to convert the .bam generated in script "01_" into something human readable. It requires a few arguments when executing from the command line:
    -i: input .bam files (see example at top of script for pointing to the right directory)
    -m: input metadata. This is the same .txt file utilized in script "01_"
    -r: reference folder. contains fasta files of reference sequences
    -n: count per million cutoff. For 1E6 reads, will disregard all sequences read less than n/1E6 times. This scales with the total number of reads mapped within the .bam file.

The script can be run multiple times with different -n values and the output will be generated separately from previous runs with the cpm value noted in the output file name.

Upon completion, script "02_" creates a new subfolder within "outputs" titled "mutation_dfs" which will ouptut a .csv with the following components for all full length sequences identified in each .bam:
    1. DNA sequence: full DNA sequence mapped to refernce
    2. Amino Acid Sequence: translation of DNA sequence in (1.)
    3. Read Count: number of times that sequence is observed in the data
    4. Normalized Read Count
    5. DNA Mutations
    6. Amino Acid Mutations
    7. Number of DNA Mutations
    8. Number of Amino Acid Mutations

Note that the DNA and Amino Acid Sequences populated will only match the "reference_size" value set in the "metadata.txt" file. This file is convenient to quickly glance at the level of diversity within your data set. This will also be the starting file for additional downstream scripts used to graph the data within.

Script "03_" will be added in a subsequent push
