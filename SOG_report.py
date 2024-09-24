import sys
import pandas as pd
from Bio import SeqIO
from Bio import SearchIO
import os
import glob
import multiprocessing as mp
from ete3 import Tree
from functools import partial
import subprocess
from pathlib import Path
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Levenshtein import distance 


new_site = '/g/data/fa63/zl1602/software/SOG/lib/python3.12/site-packages'
sys.path.insert(0, new_site)

#Here set the chromosome label mapping:
data_string = """AU3 JM1 MF1 GM1
CHR01 CHR02 CHR03 CHR01
CHR02 CHR01 CHR02 CHR02
CHR03 CHR03 CHR04 CHR03
CHR04 CHR04 CHR05 CHR04
CHR09 CHR09 CHR09 CHR05
CHR05 CHR05 CHR06 CHR06
CHR08 CHR10 CHR01 CHR07
CHR06 CHR06 CHR07 CHR08
CHR11 CHR11 CHR10 CHR09
CHR10 CHR12 CHR11 CHR10
CHR12 CHR07 CHR08 CHR11
CHR15 CHR15 CHR13 CHR12
CHR07 CHR08 CHR17 CHR13
CHR13 CHR13 CHR14 CHR14
CHR14 CHR14 CHR12 CHR15
CHR16 CHR16 CHR16 CHR16
CHR17 CHR17 CHR15 CHR17
CHR18 CHR18 CHR18 CHR18"""

# Split the string into lines
lines = data_string.strip().split('\n')

# The first line contains the column names
columns = lines[0].split()

# The rest of the lines contain the data
data = [line.split() for line in lines[1:]]

# Create the DataFrame
chr_mapping_df = pd.DataFrame(data, columns=columns)
chr_mapping_df['AU3'] = 'AU3_' + chr_mapping_df['AU3']
chr_mapping_df['JM1'] = 'JM1_' + chr_mapping_df['JM1']
chr_mapping_df['MF1'] = 'MF1_' + chr_mapping_df['MF1']
chr_mapping_df['GM1'] = 'GM1_' + chr_mapping_df['GM1']
#Convert each row as set
chr_mapping_df = chr_mapping_df.apply(lambda x: set(x), axis=1)


def verify_primer(OG_name):
    final_table = os.path.join(primersearchdir_tmp, OG_name + '.filtered.primersearch')
    if os.path.exists(final_table):
        primer_result = pd.read_csv(final_table, sep='\t', header=0)
        #Make genome column as the hitseq split by second underscore, get the first and second element
        primer_result['genome'] = primer_result['hitseq'].str.split('_').str[0] + '_' + primer_result['hitseq'].str.split('_').str[1]
        primer_count = primer_result.groupby(['primer_name', 'genome']).count().reset_index()
        #If sum of primer count is equal to the number of genomes, and every primer has exact 1 hit in each genome, then the primer is good
        primer_can = primer_count['primer_name'].unique()
        primer_file = os.path.join(workdir, 'primer_tmp', OG_name + '.primer')
        primer_file = pd.read_csv(primer_file, sep='\t', header=None)
        primer_file.columns = ['primer_name', 'fwd_seq', 'rev_seq']
        for primer in primer_can:
            QC_primer(primer_result, primer, primer_file)
    else:
        return
    return



def QC_primer(primer_result_tmp, primer_name, primer_file0):
    primer_subdf = primer_result_tmp[primer_result_tmp['primer_name'] == primer_name]
    primer_report = primer_file0[primer_file0['primer_name'] == primer_name]
    #check whether the primer has hit in all genomes, it might have multiple hits in the same genome
    if len(primer_subdf['genome'].unique()) == len(fasta_file):
        #Check whether the primer has only 1 hit in each genome
        if primer_subdf['genome'].value_counts().max() == 1:
            #find the sequence from genome
            fna_search = fna_list.copy()
            target_fna_list = {}
            for i in range(len(primer_subdf)):
                genome = primer_subdf.iloc[i, 1]
                seq = [x for x in fna_search if x.id == genome][0]
                strand = primer_subdf.iloc[i, 7]
                header = seq.id + '_' + primer_name
                start = primer_subdf.iloc[i, 2]
                end = primer_subdf.iloc[i, 4]
                if strand == '+':
                    seq = seq.seq[start:end]
                else:
                    seq = seq.seq[start:end].reverse_complement()
                target_fna_list[header] = seq
            #Now calculate the pairwise Levendstein distance
            primer_distance = {}
            keys = list(target_fna_list.keys()) 
            for i, key1 in enumerate(keys):
                for key2 in keys[i+1:]:  # Start from the next element to avoid duplicate combinations
                    distance_tmp = distance(str(target_fna_list[key1]), str(target_fna_list[key2]))
                    primer_distance[key1 + '_' + key2] = distance_tmp
            #Count the number of 0 from the primer_distance
            zero_count = sum([1 for x in primer_distance.values() if x == 0])
            max_distance = max(primer_distance.values())
            distance_df = pd.DataFrame(index=keys, columns=keys)
            for i, key1 in enumerate(keys):
                for key2 in keys[i+1:]:
                    distance_df.loc[key1, key2] = primer_distance[key1 + '_' + key2]
                    distance_df.loc[key2, key1] = primer_distance[key1 + '_' + key2]
            #Now confirming whether the target are on the same chromosome
            target_match_list = [x.split('_')[0] + '_' + x.split('_')[2] for x in keys]
            target_match_list = [x.replace('AB', '') for x in target_match_list]
            target_match_list = [x.replace('ab', '') for x in target_match_list]
            target_set = set(target_match_list)
            primer_report = primer_report.copy()
            primer_report.loc[:,'max_differences'] = max_distance
            primer_report.loc[:,'numbers_of_identical_genome'] = zero_count
            primer_report.loc[:,'min_amplicon_length'] = primer_subdf['amplicon_len'].min()
            primer_report.loc[:,'max_amplicon_length'] = primer_subdf['amplicon_len'].max()
            primer_report.loc[:,'max_missmatch_fwd'] = primer_subdf['fwd_mismatch'].max()
            primer_report.loc[:,'max_missmatch_rev'] = primer_subdf['rev_mismatch'].max()
            if len(target_set) > 4:
                primer_report.loc[:,'translocation'] = 'Yes'
            else:
                #Compare the set to chr_mapping_df
                for i in range(len(chr_mapping_df)):
                    if target_set == chr_mapping_df.iloc[i]:
                        primer_report.loc[:,'translocation'] = 'No'
                        #use the element contain 'AU3' as the chromosome name
                        chr_name = [x for x in target_set if 'AU3' in x][0]
                        primer_report.loc[:,'chromosome'] = chr_name
                        break
                    else:
                        primer_report.loc[:,'translocation'] = 'Yes'
                        primer_report.loc[:, 'chromosome'] = 'Uncertain'
            #Now write the primer_report to a file
            OG_name0 = primer_name.split('_')[0]
            primer_report.to_csv(os.path.join(primer_report_dir, OG_name0 + '.report.csv'), sep='\t', index=False, header=False, mode='a')
            distance_df.to_csv(os.path.join(primer_report_dir, primer_name + '.distance.csv'), sep='\t', index=True, header=True)
        else:
            return
    else:
        return
    

def process_files_parallel(OG_name, num_threads):
    pool = mp.Pool(num_threads)
    pool.map(verify_primer, OG_name)
    pool.close()
    pool.join()


if __name__ == '__main__':
    num_threads = 32
    workdir = os.getenv('PBS_JOBFS')
    workdir = os.path.join(workdir, 'SOG')
    fasta_folder = os.path.join(workdir, 'genome')
    fasta_file = glob.glob(fasta_folder + '/*.fa')
    primer_folder = os.path.join(workdir, 'primer_tmp')
    primersearchdir_tmp = os.path.join(workdir, 'primersearch_tmp')
    primer_candidate = glob.glob(primersearchdir_tmp + '/*.primersearch')
    OG_name = [os.path.basename(x).split('.')[0] for x in primer_candidate]
    primer_report_dir = os.path.join(workdir, 'primer_report')
    if not os.path.exists(primer_report_dir):
        os.makedirs(primer_report_dir)
    
    fna_list = []

    for fasta in fasta_file:
        fna = list(SeqIO.parse(fasta, 'fasta'))
        fna_list.extend(fna)

    process_files_parallel(OG_name, num_threads)