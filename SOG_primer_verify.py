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
emboss_path = '/g/data/fa63/zl1602/software/EMBOSS-6.6.0/emboss'
os.environ['PATH'] = emboss_path + ':' + os.environ['PATH']

try:
    result = subprocess.run(['primersearch', '-help'], check=False, capture_output=True)
    print("primersearch is available. ")
except subprocess.CalledProcessError as e:
    print("primersearch failed to run or is not available")

workdir = os.getenv('PBS_JOBFS')
workdir = os.path.join(workdir, 'SOG')
primersearchdir = os.path.join(workdir, 'primersearch')
primer_report_dir = os.path.join(workdir, 'primer_report')
if not os.path.exists(primer_report_dir):
    os.makedirs(primer_report_dir)

if not os.path.exists(primersearchdir):
    os.makedirs(primersearchdir)

fasta_folder = os.path.join(workdir, 'genome')

fasta_file = glob.glob(fasta_folder + '/*.fa')


def read_primersearch(primersearch_file, primer_file0):
    primer_result = pd.DataFrame(columns=['primer_name', 'hitseq', 'fwd_start', 'fwd_mismatch', 'rev_end', 'rev_mismatch', 'amplicon_len', 'strand'])
    with open(primersearch_file, 'r') as f:
        lines = f.readlines()
    print(f'Processing {primersearch_file}')
    #If the primersearch result is empty, then return empty dataframe
    if len(lines) > 0:
        primer = pd.read_csv(primer_file0, sep='\t', header=None)
        primer.columns = ['primer_name', 'fwd_seq', 'rev_seq']
        primer_fwd_seq = primer['fwd_seq'].tolist()
        for i in range(len(lines)):
            if lines[i].startswith('Primer name'):
                primer_name = lines[i].split()[2]
                hitseq = lines[i+2].split()[1]
                fwd_start = lines[i+4].split()[5]
                fwd_mismatch = lines[i+4].split()[7]
                f_seq = lines[i+4].split()[0]
                rev_mismatch = lines[i+5].split()[7]
                amplicon_len = lines[i+6].split()[2]
                amplicon_len = int(amplicon_len)
                if f_seq in primer_fwd_seq: 
                    rev_end = int(fwd_start) + int(amplicon_len) - 1   
                    fwd_start = int(fwd_start) -1
                    strand = '+'
                else:
                    rev_end = int(fwd_start) + int(amplicon_len) - 1
                    fwd_start = int(fwd_start) -1
                    strand = '-'
                tmp_df = pd.DataFrame({'primer_name':primer_name, 'hitseq':hitseq, 'fwd_start':fwd_start, 'fwd_mismatch':fwd_mismatch, 'rev_end':rev_end, 'rev_mismatch':rev_mismatch, 'amplicon_len':amplicon_len, 'strand':strand }, index=[0])
                primer_result = pd.concat([primer_result, tmp_df], ignore_index=True)
    return primer_result



def verify_primer(OG_name):
    #Check if the primersearch result exists
    final_table = os.path.join(primersearchdir_tmp, OG_name + '.filtered.primersearch')
    #If the primersearch result exists, then skip the primersearch
    if os.path.exists(final_table):
        print(f'{OG_name} previously done. ')
        return
    else:
        primer_file = os.path.join(workdir, 'primer_tmp', OG_name + '.primer')
        for genome in fasta_file:
            prefix = os.path.basename(genome).split('.')[0]
            prefix = prefix + '_' + OG_name
            outfile_name = os.path.join(primersearchdir, f'{prefix}.primersearch')
            input_name = f'{genome}'
            command = ['primersearch' , '-infile', primer_file, '-seqall', input_name, '-mismatchpercent', '10',  '-outfile', outfile_name]
            #print(f'Running primersearch for {prefix}')
            subprocess.run(command, check=True, capture_output=True)
        #Now checking the primersearch result
        primersearch_file = glob.glob(f'{primersearchdir}/*_{OG_name}.primersearch')
        for i in range(len(primersearch_file)):
            if i == 0:
                primer_result = read_primersearch(primersearch_file[i], primer_file)
            else:
                primer_tmp = read_primersearch(primersearch_file[i], primer_file)
                primer_result = pd.concat([primer_result, primer_tmp], ignore_index=True)
        for file_path in primersearch_file:
            os.remove(file_path)
        #Filter the primer result
        primer_result = primer_result[primer_result['amplicon_len'] <= 5000]
        #Remove if the hitseq contains 'CTG' or 'SCF'
        primer_result = primer_result[~primer_result['hitseq'].str.contains('CTG')]
        primer_result = primer_result[~primer_result['hitseq'].str.contains('SCF')]
        primer_result['genome'] = primer_result['hitseq'].str.replace('_CHR*', '')
        #Group by primer_name and hitseq
        primer_result.to_csv(final_table, sep='\t', index=False)
        #copy the primer.tsv to the g/data
        command = ['cp', final_table, '/g/data/fa63/zl1602/primersearch_tmp/']
        subprocess.run(command, capture_output=True)
        print(f'{OG_name} is copied to destination. ')
  #  fna0_list = fna_list.copy()
    return
    

def process_files_parallel(OG_name, num_threads):
    pool = mp.Pool(num_threads)
    pool.map(verify_primer, OG_name)
    pool.close()
    pool.join()


if __name__ == '__main__':
    num_threads = 128
    primer_folder = os.path.join(workdir, 'primer_tmp')
    primer_candidate = glob.glob(primer_folder + '/*.primer')
    OG_name = [os.path.basename(x).split('.')[0] for x in primer_candidate]
    primer_verify_dir = os.path.join(workdir, 'primer_verify')
    if not os.path.exists(primer_verify_dir):
        os.makedirs(primer_verify_dir)

    primersearchdir_tmp = os.path.join(workdir, 'primersearch_tmp')
    if not os.path.exists(primersearchdir_tmp):
        os.makedirs(primersearchdir_tmp)
    process_files_parallel(OG_name, num_threads)
    
    