#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on July 17 2025
@author: jhambale
example run: python 03_manzanaa_plots.py \
-i ../../manzanaa_data/manzanaa_outputs/fl_001/mutation_dfs/*mutation_analysis.csv \
-m ../../manzanaa_data/ngs_raw/fl_001/demultiplex/fl_001_demux_metadata.txt \
-r ../../manzanaa_data/ngs_raw/fl_001/references/ \
-c 25 -s 40 -t 10
"""
# activate virtual enviroment before running script
# source activate minibinders

import pandas as pd
import pysam
from Bio import SeqIO
import glob
import numpy as np
from itertools import combinations
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

from utils import *
import argparse

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', nargs='+', help='List of fastq files to parse')
    parser.add_argument('-m', help='Path to metadata')
    parser.add_argument('-r', help='Path to reference database')
    parser.add_argument('-c', help='Select CPM of interest from script 02')
    parser.add_argument('-s', help='Desired max scale of mutations')
    parser.add_argument('-t', help='Desired number of ticks per graph')

    args = parser.parse_args()


    # Set the Seaborn style and context
    sns.set(style="white")

    # Set global font size and family
    plt.rcParams.update({
        'font.size': 16,
        'font.family': 'Arial'
    })
    # Set the Seaborn style
    sns.set_style("ticks")


    mutation_df_all = args.i #glob.glob('../../manzanaa_data/manzanaa_outputs/maa_001/alignments/*.bam') # input
    metadata_path = args.m #'../../manzanaa_data/ngs_raw/maa_001/demultiplex/maa_001_demux_metadata.txt' # input
    ref_dir = args.r #'../../manzanaa_data/ngs_raw/maa_001/references/' # input
    cpm = int(args.c)
    xupperlim = float(args.s) #static number of desired DNA mutation (e.g 20)
    tick_num = float(args.t) #desired number of ticks per graph

    mutation_df_paths = [path for path in mutation_df_all if f'cpm{cpm}' in path]
    
    # mutation_df_paths = glob.glob('../../manzanaa_data/manzanaa_outputs/maa_001/mutation_dfs/*mutation_analysis.csv')
    # metadata_path = '../../manzanaa_data/ngs_raw/maa_001/demultiplex/maa_001_demux_metadata.txt' # input
    metadata_df = metadata_df = pd.read_csv(metadata_path,
                                      sep='\t', index_col=None)
    # ref_dir ='../../manzanaa_data/ngs_raw/maa_001/references/' # input

    for mutation_path in mutation_df_paths:

        tick_scale = xupperlim/tick_num
        aa_xupperlim = xupperlim


    #     mutation_df_path = mutation_path
        mutation_df = pd.read_csv(mutation_path, index_col=0)



        sample_name = mutation_path.split('/')[-1].replace(f'cpm{cpm}_mutation_analysis.csv','')
        outdir = '/'.join(mutation_path.replace('mutation_dfs', 'mutation_figs').split('/')[:-1]) + '/'
        make_dir(outdir)
        outdir_sample = outdir + sample_name +'/'
        make_dir(outdir_sample)

        outname = outdir_sample + sample_name

        metadata_samp = metadata_df.copy()[metadata_df['sample_id'].str.contains(sample_name)]
        reference_name = list(metadata_samp['reference'])[0]

        num_sequences = len(mutation_df)

        if '.fasta' not in reference_name:
            reference_name = reference_name + '.fasta'
        ref_path = ref_dir + reference_name

        reference_size = metadata_samp['reference_size']#int(metadata_samp['reference_size'])
        with open(ref_path, "rt") as reference_reads:
            reference = SeqIO.parse(reference_reads, "fasta")
            for reads in reference:
                reference_sequence_dna = reads.seq
        reference_sequence_aa = str(reference_sequence_dna.translate())

        mutation_clustered_df = mutation_df.copy().drop_duplicates(subset='aa_sequence').reset_index(drop=True)
        # print(len(mutation_df))
        # print(len(mutation_clustered_df))
        dna_hamming = list(mutation_df['number_dna_mutations'])
        aa_hamming = list(mutation_clustered_df['number_aa_mutations'])

        dna_sequences = list(mutation_df['dna_sequence'])
        aa_sequences = list(mutation_clustered_df['aa_sequence'])

        hamming_pairwise_dna = []
        hamming_pairwise_aa = []
        for dna_seq_combo in combinations(dna_sequences, 2):
            hamming_pairwise_dna.append(hamming_distance(dna_seq_combo[0], dna_seq_combo[1]))
        for aa_seq_combo in combinations(aa_sequences, 2):
            hamming_pairwise_aa.append(hamming_distance(aa_seq_combo[0], aa_seq_combo[1]))
       

        g = sns.histplot(data=mutation_df, x="number_dna_mutations",
                    binwidth=1,stat='probability',
                     edgecolor='black', linewidth=0.5)
        ax = plt.gca()
        g.xaxis.set_major_locator(MaxNLocator(integer=True))


        ticks = ax.get_xticks()
        print(ticks)
        print(ticks[-1])
        # set xmax, linspace from 0, add offset as percentage of range
        if ticks[-1] > xupperlim:
            scaled_ticks = np.linspace(-(tick_scale),xupperlim, num=int(((xupperlim/tick_scale)+2)),endpoint=True)
            print(scaled_ticks)
            new_ticks = scaled_ticks[:-1] + 0.5
        else:
            new_ticks = ticks[:-1] + 0.5

        
        print(new_ticks)
        print(new_ticks[0],new_ticks[-1])
        #reset new ticks
        g.set_xticks(new_ticks)
        g.set_xticklabels([str(int(t)) for t in new_ticks])

        g.set_xlabel('dna mutations')
        g.set_ylabel(f'frequency (n={num_sequences} sequences)')
        g.spines['right'].set_visible(False)
        g.spines['top'].set_visible(False)
        for spine in ['left', 'bottom']:
            g.spines[spine].set_linewidth(0.5)
        g.tick_params(width=0.5)
        g.set_xlim([new_ticks[0]+0.5, new_ticks[-1]-0.5])


        # Show/save the plot
        plt.savefig(f'{outname}cpm{cpm}_dna_dist_wt.png', dpi=400)
        # plt.savefig(f'{outname}_dna_dist_wt.pdf', dpi=400)
        # plt.show()
        plt.close()

        g = sns.histplot(data=mutation_df, x="number_aa_mutations",
                     binwidth=1,stat='probability',
                     edgecolor='black', linewidth=0.5)
        ax = plt.gca()
        # Adjust the ticks
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))


        ticks = ax.get_xticks()
        print(ticks)
        print(ticks[-1])
        #set xmax, linspace from 0, add offset as percentage of range
        if ticks[-1] > aa_xupperlim:
            scaled_ticks = np.linspace(-(tick_scale),aa_xupperlim, num=int(((aa_xupperlim/tick_scale)+2)),endpoint=True)
            print(scaled_ticks)
            new_ticks = scaled_ticks[:-1] + 0.5
        else:
            new_ticks = ticks[:-1] + 0.5
            
        g.set_xticks(new_ticks)
        g.set_xticklabels([str(int(t)) for t in new_ticks])

        g.set_xlabel('amino acid mutations')
        g.set_ylabel(f'frequency (n={num_sequences} sequences)')
        g.spines['right'].set_visible(False)
        g.spines['top'].set_visible(False)
        for spine in ['left', 'bottom']:
            g.spines[spine].set_linewidth(0.5)
        g.tick_params(width=0.5)
        g.set_xlim([new_ticks[0]+0.5, new_ticks[-1]-0.5])

        # Show the plot
        plt.savefig(f'{outname}cpm{cpm}_aa_dist_wt.png', dpi=400)
        # plt.savefig(f'{outname}_aa_dist_wt.pdf', dpi=400)
        plt.close()
        
        # Create a seaborn histogram
        g = sns.histplot(data=hamming_pairwise_dna,
                         binwidth=1,stat='probability',
                         edgecolor='black', linewidth=0.5)
        ax = plt.gca()
        # Adjust the ticks
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))


        ticks = ax.get_xticks()
        print(ticks)
        print(ticks[-1])
        # set xmax, linspace from 0, add offset as percentage of range
        if ticks[-1] > xupperlim:
            scaled_ticks = np.linspace(-(tick_scale),xupperlim, num=int(((xupperlim/tick_scale)+2)),endpoint=True)
            print(scaled_ticks)
            new_ticks = scaled_ticks[:-1] + 0.5
        else:
            new_ticks = ticks[:-1] + 0.5


        g.set_xticks(new_ticks)
        g.set_xticklabels([str(int(t)) for t in new_ticks])

        g.set_xlabel('dna mutations')
        g.set_ylabel(f'frequency (n={num_sequences} sequences)')
        g.spines['right'].set_visible(False)
        g.spines['top'].set_visible(False)
        for spine in ['left', 'bottom']:
            g.spines[spine].set_linewidth(0.5)
        g.tick_params(width=0.5)
        g.set_xlim([new_ticks[0]+0.5, new_ticks[-1]-0.5])

        # Show the plot
        plt.savefig(f'{outname}cpm{cpm}_dna_dist_pairwise.png', dpi=400)
        # plt.savefig(f'{outname}_dna_dist_pairwise.pdf', dpi=400)

        plt.close()

        # Create a seaborn histogram
        g = sns.histplot(data=hamming_pairwise_aa,
                         binwidth=1,stat='probability',
                         edgecolor='black', linewidth=0.5)
        ax = plt.gca()
        # Adjust the ticks
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))


        ticks = ax.get_xticks()
        print(ticks)
        print(ticks[-1])
        #set xmax, linspace from 0, add offset as percentage of range
        if ticks[-1] > xupperlim:
            scaled_ticks = np.linspace(-(tick_scale),aa_xupperlim, num=int(((aa_xupperlim/tick_scale)+2)),endpoint=True)
            print(scaled_ticks)
            new_ticks = scaled_ticks[:-1] + 0.5
        else:
            new_ticks = ticks[:-1] + 0.5
        
        
        g.set_xticks(new_ticks)
        g.set_xticklabels([str(int(t)) for t in new_ticks])

        g.set_xlabel('amino acid mutations')
        g.set_ylabel(f'frequency (n={num_sequences} sequences)')
        g.spines['right'].set_visible(False)
        g.spines['top'].set_visible(False)
        for spine in ['left', 'bottom']:
            g.spines[spine].set_linewidth(0.5)
        g.tick_params(width=0.5)
        g.set_xlim([new_ticks[0]+0.5, new_ticks[-1]-0.5])


        # Show the plot
        plt.savefig(f'{outname}cpm{cpm}_aa_dist_pairwise.png', dpi=400)
        # plt.savefig(f'{outname}_aa_dist_pairwise.pdf', dpi=400)
        plt.close()

        # Create a seaborn histogram
        g = sns.histplot(data=hamming_pairwise_aa,
                         binwidth=1,stat='probability',
                         edgecolor='black', linewidth=0.5)
        sns.histplot(data=mutation_df, x="number_aa_mutations",
                         binwidth=1,stat='probability',
                         edgecolor='black', linewidth=0.5, color='gray')
        ax = plt.gca()
        # Adjust the ticks
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))



        ticks = ax.get_xticks()
        print(ticks)
        print(ticks[-1])
        #set xmax, linspace from 0, add offset as percentage of range
        if ticks[-1] > xupperlim:
            scaled_ticks = np.linspace(-(tick_scale),aa_xupperlim, num=int(((aa_xupperlim/tick_scale)+2)),endpoint=True)
            print(scaled_ticks)
            new_ticks = scaled_ticks[:-1] + 0.5
        else:
            new_ticks = ticks[:-1] + 0.5
        
        
        # ticks = g.get_xticks()
        # new_ticks = ticks[:-1] + 0.5
        g.set_xticks(new_ticks)
        g.set_xticklabels([str(int(t)) for t in new_ticks])

        g.set_xlabel('amino acid mutations')
        g.set_ylabel(f'frequency (n={len(mutation_clustered_df)} sequences)')
        g.spines['right'].set_visible(False)
        g.spines['top'].set_visible(False)
        for spine in ['left', 'bottom']:
            g.spines[spine].set_linewidth(0.5)
        g.tick_params(width=0.5)
        g.set_xlim([new_ticks[0]+0.5, new_ticks[-1]-0.5])


        # Show the plot
        plt.savefig(f'{outname}cpm{cpm}_ham_pw.png', dpi=400)
        # plt.savefig(f'{outname}_ham_pw.pdf', dpi=400)
        plt.close()
        print(np.mean(hamming_pairwise_aa))
        print(np.mean(aa_hamming))

if __name__ == "__main__":
    main()
