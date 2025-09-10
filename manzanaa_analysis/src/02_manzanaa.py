#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on July 02 2025
@author: jhambale
example run 1: python 02_manzanaa.py  -i ../../manzanaa_data/manzanaa_outputs/fl_001/alignments/*.bam \
-m ../../manzanaa_data/ngs_raw/fl_001/demultiplex/fl_001_demux_metadata.txt \
-r ../../manzanaa_data/ngs_raw/fl_001/references/ \
-n 5
example run 2: python 02_mutation_analysis.py  -i ../../manzanaa_data/manzanaa_outputs/fl_002/alignments/*.bam -m ../../manzanaa_data/ngs_raw/fl_002/demultiplex/fl_002_demux_metadata.txt -r ../../manzanaa_data/ngs_raw/fl_002/references/ -n 50

"""
# import relevant packages

import pandas as pd
import pysam
from collections import Counter
from Bio import SeqIO
import glob

from utils import *

import argparse

def extract_sequence(read, limits):

    """
    Given a pysam read and the index limits of a desired region, extract out the
    sequence of the desired region of the read, adjusting for any insertions or
    deletions in the alignment. written by V. Hu. Adapted by M. Alcantar
    to work with versions below python 3.10 (i.e., removed the match operator)

    PARAMETERS
    -----------
    read: pysam AlignedSegment object
        basically a record from a .bam file
    limits: 2-d list ([start, end])
        a list of two ints denoting index limits for the extraction

    RETURNS
    --------
        a string denoting the sequence extracted from the reads
    """
    index = 0  # Pointer to position on the query
    q_start, q_end = limits  # start and end positions on the reference
    q_start -= read.reference_start
    q_end -= read.reference_start
    for operation, positions in read.cigartuples:
        if operation in [0, 7, 8]:  # Consume both Q and R
            index += positions
        elif operation == 1:  # Consume only Q
            # move bookmarks up to account for x extra query positons
            q_start += (index <= q_start) * positions
            q_end += (index < q_end) * positions
            index += positions
        elif operation in [2, 3]:  # Consume only R
            # move bookmarks down to account for x ref skips
            q_end -= (index < q_end) * min(positions, q_end - index)
            q_start -= (index < q_start) * min(positions, q_start - index)
        elif operation == 4:  # Consume only Q, but only at start and end
            q_start += (index <= q_start or q_start <= 0) * positions
            q_start += (q_start <= 0 or (index < q_start and index + positions > q_start)) * (
                        index + positions - q_start)
            q_end += (index <= q_end or q_end <= 0) * positions
            index += positions
            q_end += (index - positions < q_end and index > q_start) * (index - q_end - positions)
    return(read.seq[max(q_start, 0): q_end])

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', help='List of alignment files to parse')
    parser.add_argument('-m', help='Path to metadata')
    parser.add_argument('-r', help='Path to reference database')
    parser.add_argument('-n', default=50, help='counts per million')

    args = parser.parse_args()

    alignment_files = args.i #glob.glob('../../minibinders_data/minibinders_outputs/fl_001/alignments/*.bam') # input
    metadata_path = args.m #'../../minibinders_data/ngs_raw/fl_001/demultiplex/fl_001_demux_metadata.txt' # input
    ref_dir = args.r#'../../minibinders_data/ngs_raw/fl_001/references/' # input
    num_seq_filter = int(args.n)#100 #input

    # in case you forget to add the "/" at the end of your references input ;)
    if ref_dir[-1] != '/':
        ref_dir = ref_dir + '/'

    metadata_df = pd.read_csv(metadata_path,
                              sep='\t', index_col=None)

    for alignment in alignment_files:

        print('----------------------------')
        print(f'parsing: {alignment}')
        print('----------------------------')

        sample_name = alignment.split('/')[-1].replace('.bam','')
        metadata_samp = metadata_df.copy()[metadata_df['sample_id'].str.contains(sample_name)]
        reference_name = list(metadata_samp['reference'])[0]

        if '.fasta' not in reference_name:
            reference_name = reference_name + '.fasta'
        ref_path = ref_dir + reference_name

        print(metadata_samp['reference_size'])
        reference_size = int(metadata_samp['reference_size'])
        samfile = pysam.AlignmentFile(alignment, "rb")
        with open(ref_path, "rt") as reference_reads:
            reference = SeqIO.parse(reference_reads, "fasta")
            for reads in reference:
                reads = reads.upper()
                reference_sequence_dna = reads.seq
        reference_sequence_aa = str(reference_sequence_dna.translate())

        # Open the BAM file
        sequences = []
        aa_sequences = []
        print('extracting sequences')
        for read in samfile:
            # extract sequences and only take sequences if they match the expected size
            try:
                sequence = extract_sequence(read, [0,reference_size])
                if len(sequence) == reference_size and 'N' not in sequence:
                    sequences.append(sequence)
            except:
                pass

        # filter sequences:
        # find sequences that appear cpm times
        print('filtering sequences')
        sequences_counter = Counter(sequences)
        # print(sequences_counter)
        print(f"{sum(sequences_counter.values())} total sequences for {sample_name}")

        mutation_analysis_dict_list = []
        num_sequences_filter = int(num_seq_filter/10**6 *sum(sequences_counter.values())) #counts per million
        print(f"read count threshold: {num_sequences_filter}")
        
        filtered_sequences = {seq: count for seq, count in sequences_counter.items() if count >= num_sequences_filter}

        filtered_sequences_dict = filtered_sequences
        seq_count_total = sum(filtered_sequences_dict.values())
        print('finding mutations')
        for dna_variant in filtered_sequences_dict:
            sequence_dna_tmp = dna_variant
            sequence_aa_tmp = translate_sequence(dna_variant)
            sequence_cnt_tmp = filtered_sequences_dict[dna_variant]
            normalized_read_count = sequence_cnt_tmp/seq_count_total

            dna_muts = find_mutations(reference_sequence_dna, sequence_dna_tmp)
            number_dna_mutations = len(dna_muts)
            aa_muts = find_mutations(reference_sequence_aa, sequence_aa_tmp)
            number_aa_mutations = len(aa_muts)
            mutation_analysis_dict_list.append({'dna_sequence':sequence_dna_tmp,
                                               'aa_sequence': sequence_aa_tmp,
                                               'read_count': sequence_cnt_tmp,
                                               'normalized_read_count': normalized_read_count,
                                               'dna_mutations': dna_muts,
                                               'aa_mutations': aa_muts,
                                               'number_dna_mutations':number_dna_mutations,
                                               'number_aa_mutations': number_aa_mutations})

        mutations_df_outdir = '/'.join(alignment.split('/')[:-1]).replace('alignments', 'mutation_dfs/')
        make_dir(mutations_df_outdir)
        mutations_df_outname = sample_name + f'cpm{str(num_seq_filter)}_mutation_analysis.csv'
        mutations_df_outpath =mutations_df_outdir +  mutations_df_outname
        mutations_df = pd.DataFrame(mutation_analysis_dict_list)
        mutations_df = mutations_df.sort_values(by='read_count', ascending=False).reset_index(drop=True)
        mutations_df.to_csv(mutations_df_outpath)

if __name__ == "__main__":
    main()
