#!/bin/bash

######################################################

# This script performs the following:
# i) FastQC - extracts quality metrics from raw fastq files
# ii) Fastp - performs quality filtering and trimming on raw reads
# iii) Seqkit - extracts general summary statistics for sequencing files
# iv) minimap2 - performs alignment to a reference sequence
# outputs are written to ../../manzanaa_data/manzanaa_outputs directory
#
# this script is best used with single read data
# the filtering metrics are loosened to account for Premium PCR's lower quality
#
### PARAMETERS
## $1: metadata in .txt format (tab separated) see example
# headers should include:
# sample_number (e.g. 1); experiment_name (e.g. maa-001);
# sample_name (e.g. sample1_R1_001.fastq.gz) - should match .fastq.gz file names;
# barcode_f (e.g., atgc); barcode_r (e.g., atgc) if applicable
## $2: experiment_name given explicity as a separate argument. should match the experiment name given in metadata (used for md5sum check)
#
## $3+ list of fastq or fastq.gz files to be run through (loaded as wildcard)
# Only specify one of each read pair; the other
# filename is assumed
# Example 1: bash 01b_filter_nanopore_reads.sh
# ../../manzanaa_data/ngs_raw/fl_001/demultiplex/fl_001_demux_metadata.txt fl_001
# ../../manzanaa_data/ngs_raw/jh_001/demultiplex/*.fastq.gz
#
### OUTPUTS
## FastQC
## ${fastqc_out_path}.html : fastqc report -- web link
## ${fastqc_out_path}.zip : fastqc report -- detailed files
#
## Fastp
## ${fastp_out_path}_out.fastq.gz : quality filtered read (lenient to account for nanopore accuracy)
#
## Seqkit
## ${seqkit_stats_out_path}_filtering_stats.txt: statistics on number of reads
# containing barcode
#
## minimap2
## ${alignment_out_path}.bam: alignment file
#
## NOTES
# this script assumes the file name ends in R[12]_001.fastq.gz
# any string is allowed before R[12], as long as the forward and reverse read
# prefix is the same

######################################################

METADATA=$1
experiment_name=$2
shift 2
# shift

echo $experiment_name

# check for appropriate download of fastq files through md5check
md5sum $@ > ../../manzanaa_data/ngs_raw/$experiment_name/demultiplex/md5sum_check.txt
if ! md5sum -c ../../manzanaa_data/ngs_raw/$experiment_name/demultiplex/md5sum_check.txt > ../../manzanaa_data/ngs_raw/$experiment_name/demultiplex/md5_check.log 2>&1; then
    echo "Error: corrupted file download."
    cat ../../manzanaa_data/ngs_raw/$experiment_name/demultiplex/md5_check.log
    exit 1
fi



# check to make sure files have the expected formats
for fastq in $@
do
	if [ -z $(basename $fastq | grep -i .fastq) ]
	then
		echo $(basename $fastq) "does not have .fastq suffix - aborting"
		exit 1
	fi

done
echo "all files are the correct format. continuing..."

# loop through all fastq files and extract file names and metadata required for
# writing outputs to correct folder
for fastq in "$@"
do
	fname=$(basename $fastq)
	dname=$(dirname $fastq)
	fpath=$dname/${fname%_R[12]_001.fastq*}
	experiment_name=$(cat $METADATA | grep $fname | cut -f2)
	reference_name=$(cat $METADATA | grep $fname | cut -f9)
    reference_size=$(cat $METADATA | grep $fname | cut -f10)
	flank_f=$(cat $METADATA | grep $fname | cut -f7)
	flank_r=$(cat $METADATA | grep $fname | cut -f8)

  echo Experiment name $experiment_name
    fastqc_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/fastqc_output/
    fastp_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/fastp_output/${fname%_R[12]_001.fastq*}
    # merge_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/merged_reads/${fname%_R[12]_001.fastq*}
	seqkit_stats_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/seqkit_stats/${fname%_R[12]_001.fastq*}
	seqkit_filtered_read_variants_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/merged_variants_filtered_reads/${fname%_R[12]_001.fastq*}
	alignment_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/alignments/${fname%_R[12]_001.fastq*}
  unmapped_out_path=../../manzanaa_data/manzanaa_outputs/$experiment_name/alignments/unmapped/${fname%_R[12]_001.fastq*}
	reference_path=../../manzanaa_data/ngs_raw/$experiment_name/references/${reference_name%.fasta}


	# define prefix for all output directories
	output_dir_prefix=../../manzanaa_data/manzanaa_outputs/$experiment_name
	echo $output_dir_prefix
	if [ ! -e ${merge_out_path}_out_merged.assembled.fastq ]
	then

	  # create output directories if they don't exist already
		if [[ "$fastq" == "$1" ]]
			then

				# create general output folder for current experiment
				if [ ! -d "$output_dir_prefix" ]
			  then
					echo Making directory $output_dir_prefix
			    mkdir $output_dir_prefix
				fi

				# output directory for fastqc results
				if [ ! -d "$output_dir_prefix/fastqc_output/" ]
			  then
					echo Making directory $output_dir_prefix/fastqc_output/
			    mkdir $output_dir_prefix/fastqc_output/
				fi

				# output directory for fastp results
	      if [ ! -d "$output_dir_prefix/fastp_output/" ]
			  then
					echo Making directory $output_dir_prefix/fastp_output/
			    mkdir $output_dir_prefix/fastp_output/
				fi
				# output directory for fastp results
	      if [ ! -d "$output_dir_prefix/alignments/" ]
			  then
					echo Making directory $output_dir_prefix/alignments/
			    mkdir $output_dir_prefix/alignments/
          mkdir $output_dir_prefix/alignments/unmapped/
				fi

				# # output directory for merged reads
	      # if [ ! -d "$output_dir_prefix/merged_reads/" ]
			  # then
				# 	echo Making directory $output_dir_prefix/merged_reads/
			  #   mkdir $output_dir_prefix/merged_reads/
				# fi

				# output directory for merged reads
	      if [ ! -d "$output_dir_prefix/merged_variants_filtered_reads/" ]
			  then
					echo Making directory $output_dir_prefix/merged_variants_filtered_reads/
			    mkdir $output_dir_prefix/merged_variants_filtered_reads/
				fi

				# output directory for seqkit stats
	      if [ ! -d "$output_dir_prefix/seqkit_stats/" ]
			  then
					echo Making directory $output_dir_prefix/seqkit_stats/
			    mkdir $output_dir_prefix/seqkit_stats/
				fi
			fi


    	 fi

      # run fastqc on forward and reverse reads
	  fastqc ${fpath}_R1_001.fastq.gz -o $fastqc_out_path
	  # fastqc ${fpath}_R2_001.fastq.gz -o $fastqc_out_path

	  # run fastp (with default settings) on forward and reverse reads
		# save quality report
	  fastp \
	      -i ${fpath}_R1_001.fastq.gz \
	      -o ${fastp_out_path}_out_R1.fastq \
	      -j ${fastp_out_path}.json \
	      -h ${fastp_out_path}.html \
          -q 15 \
          -u 15 \
          --length_required 50


	  # run merging (via NGmerge) on quality filtered forward and reverse reads

      # NGmerge \
      #     -1 "${fastp_out_path}_out_R1.fastq.gz" \
      #     -2 "${fastp_out_path}_out_R2.fastq.gz" \
      #     -o "${merge_out_path}_out_ngmerged.assembled.fastq" \
      #         -n 30 \
      #         -f ${merge_out_path}_out_ngmerged.unassembled.fastq \
      #     -v \
      #     2> ${merge_out_path}_out_ngmerge_stats.txt
      #
      # gunzip "${merge_out_path}_out_ngmerged.assembled.fastq"

	### FILTER FASTQ FILES USING SEQKIT
	echo Beginning filtering for $fname

	# extract reads that contain the desired barcode configuration
	# (e.g., attcNNNNNNNNNcAattcNNNNNNNNNcAATT for doubles)
	cat ${fastp_out_path}_out_R1.fastq | seqkit stats > ${seqkit_stats_out_path}_filtering_stats.txt
	cat ${fastp_out_path}_out_R1.fastq | seqkit amplicon -m 0 -P -s -F $flank_f -R $flank_r -r 21:-21 --bed > ${seqkit_filtered_read_variants_out_path}_variants.txt

    minimap2 \
    -ax sr \
    -t 7 \
    -B 2 \
    -O 2,24 \
    -k 11 \
    $reference_path.fasta \
    ${fastp_out_path}_out_R1.fastq | samtools view -bS - > ${alignment_out_path}.bam


    samtools sort -o ${unmapped_out_path}.sorted.bam ${alignment_out_path}.bam
    samtools index ${unmapped_out_path}.sorted.bam
    samtools stats ${unmapped_out_path}.sorted.bam > ${unmapped_out_path}_sortedstats.txt


    samtools view -b -f 4 ${unmapped_out_path}.sorted.bam > ${unmapped_out_path}_unmapped.bam
    samtools fastq ${unmapped_out_path}_unmapped.bam > ${unmapped_out_path}_unmapped.fastq
    fastqc ${unmapped_out_path}_unmapped.fastq -o ../../manzanaa_data/manzanaa_outputs/$experiment_name/alignments/unmapped/


			echo _____________new sample____________________________

done
