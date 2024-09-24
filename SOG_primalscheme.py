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

new_site = '/g/data/fa63/zl1602/software/SOG/lib/python3.12/site-packages'
sys.path.insert(0, new_site)
try:
    result = subprocess.run(['primalscheme', 'multiplex', '--help'], check=False, capture_output=True)
    print("primalscheme is available. ")
except subprocess.CalledProcessError as e:
    print("primalscheme failed to run or is not available")

def count_gaps_before_nucleotide(sequence):
    gap_count = 0
    for char in sequence:
        if char == '-':
            gap_count += 1
        elif char in 'ATCGU':
            break
    return gap_count

def count_gaps_after_nucleotide(sequence):
    gap_count = 0
    for char in sequence[::-1]:
        if char == '-':
            gap_count += 1
        elif char in 'ATCGU':
            break
    return gap_count

def primer_detection(OG_name):
    file_name = OG_name + '.clustalo.aln'
    #Read the alignment file
    alignment = AlignIO.read(f'{workdir}/{file_name}', 'fasta')
    #Count the number of gaps before the first nucleotide
    start_gap = max([count_gaps_before_nucleotide(x.seq) for x in alignment])
    end_gap = max([count_gaps_after_nucleotide(x.seq) for x in alignment])
    seqlen = len(alignment[0].seq)
    end_gap = seqlen - end_gap
    #Get the maximum gap count at the beginning and end
    alignment = alignment[:, start_gap:end_gap]
    seq_list = [x.seq for x in alignment]
    file_name = OG_name + '.trim.aln'
    #Now write the alignment to a new file named as trimmed alignment
    with open(f'{workdir}/{file_name}', 'w') as f:
        for i in range(len(alignment)):
            f.write(f'>{alignment[i].id}\n')
            f.write(f'{seq_list[i]}\n')
    #Print the length of sequences in the alignment
#    amplicon_force_start = 200
#    amplicon_force_end = len(alignment[0].seq) -200
    #Set the exact position of the amplicon
    #If workdir does not exist, create it
    if not os.path.exists(os.path.join(workdir, OG_name)):
        os.makedirs(os.path.join(workdir, OG_name))
#    command = ['primalscheme', 'multiplex', '-o', f'{workdir}/{OG_name}', '-a', str(amplicon_force_start), '-a', str(amplicon_force_end), '--name', OG_name, '--force', f'{workdir}/{file_name}']
    #here if we force the amplicon size
    command = ['primalscheme', 'multiplex', '-o', f'{workdir}/{OG_name}', '-a', str(1500) , '--name', OG_name, '--force', f'{workdir}/{file_name}']
    subprocess.run(command, capture_output=True)
    #Remove all files except the primer.tsv
    out_tmp = glob.glob(f'{workdir}/{OG_name}/{OG_name}.*')
    #Drop the primer.tsv from the list
    out_tmp = [x for x in out_tmp if 'primer.tsv' not in x]
    #Remove all files in the list
    for x in out_tmp:
        os.remove(x)
    #Find whether the primer.tsv is available
    if os.path.exists(f'{workdir}/{OG_name}/{OG_name}.primer.tsv'):
        primalscheme_output = f'{workdir}/{OG_name}/{OG_name}.primer.tsv'
        primer_tsv = pd.read_csv(primalscheme_output, sep='\t', header=0)
        #Now convert primer to format that can be read by primersearch
        primer_left = primer_tsv[primer_tsv['name'].str.contains('LEFT')]
        primer_right = primer_tsv[primer_tsv['name'].str.contains('RIGHT')]
        primer_left = primer_left[['name', 'seq']]
        primer_left.columns = ['name', 'left_seq']
        primer_right = primer_right[['name', 'seq']]
        primer_right.columns = ['name', 'right_seq']
        primer_left['name']= primer_left['name'].str.replace('_LEFT', '')
        primer_right['name']= primer_right['name'].str.replace('_RIGHT', '')
        primer_df = primer_left.merge(primer_right, on='name')
        primer_df.to_csv(os.path.join(primer_tmp, OG_name + '.primer'), sep='\t', index=False, header=False)
        print(f'{OG_name} is done')
    return

def process_files_parallel(OG_name, num_threads):
    pool = mp.Pool(num_threads)
    pool.map(primer_detection, OG_name)
    pool.close()
    pool.join()




if __name__ == '__main__':
    workdir = os.getenv('PBS_JOBFS')
    workdir = os.path.join(workdir, 'SOG')
    output_folder = os.path.join(workdir, 'Primer')
    fasta = os.path.join(workdir, 'genome')
    primer_tmp = os.path.join(workdir, 'primer_tmp')
    workdir = os.path.join(workdir, 'MSA')
    os.chdir(workdir)
    num_threads = 32
    WGA_list = glob.glob(workdir + '/*.clustalo.aln')
    OG_name = [os.path.basename(x).split('.')[0] for x in WGA_list]
    
    if not os.path.exists(primer_tmp):
        os.makedirs(primer_tmp)
    process_files_parallel(OG_name, num_threads)