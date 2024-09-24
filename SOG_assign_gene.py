#!/usr/bin/env python3
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

new_site = '/g/data/fa63/zl1602/software/SOG/lib/python3.12/site-packages'
sys.path.insert(0, new_site)
#Add clustalw to the path
clustalw2_path = '/g/data/fa63/zl1602/software/clustalw-2.1/src'
os.environ['PATH'] = clustalw2_path + ':' + os.environ['PATH']
#Test the clustalw
command = ['clustalw2', '-help']
#print path
try:
    # Run the command and capture the output
    result = subprocess.run(command, check=False, capture_output=True)
    print("clustalw2 is available. ")
except subprocess.CalledProcessError as e:
    print("clustalw2 failed to run or is not available")
try:
    result = subprocess.run(['clustalo', '-help'], check=False, capture_output=True)
    print("clustalo is available. ")
except subprocess.CalledProcessError as e:
    print("clustalo failed to run or is not available")

OG_tsv = pd.read_csv("OG_output/Orthogroups.tsv", sep='\t', header=0)


def find_candidate(OG_name, output_folder):
###############Import gff here
    gff_df_cds = gff_df[gff_df['type'] == 'CDS']
    gff_df_gene = gff_df[gff_df['type'] == 'gene']
    gff_df_gene.loc[:, 'attributes'] = gff_df_gene['attributes'].str.split(';').str[0].str.split('=').str[1]
    gff_df_cds.loc[:, 'attributes'] = gff_df_cds['attributes'].str.split('.').str[0].str.split('=').str[1]
    gff_df_cds = gff_df_cds.groupby(['seqid', 'source', 'type', 'strand','attributes']).agg({'start':'min', 'end':'max'}).reset_index()
###############Import fasta here
    fna0_list = fna_list
###############Import OG here
    OG_tmp = OG_tsv[OG_tsv['Orthogroup'] == OG_name]
    cds_names = OG_tmp.iloc[0, 1:].tolist()
    gff_df_cds = gff_df_cds[gff_df_cds['attributes'].isin(cds_names)]
    gene_names = [x.split('-')[0] for x in cds_names]
    gff_df_gene = gff_df_gene[gff_df_gene['attributes'].isin(gene_names)]
    #Now get the sequence
    OG_seq = {}
    for i in range(len(gff_df_cds)):
        line = gff_df_cds.iloc[i]
        seqid = line['seqid']
        start = line['start']
        end = line['end']
        strand = line['strand']
        geneid = line['attributes']            
        seq0 = [x for x in fna0_list if x.id == seqid][0]
        if strand == '+':
            start = start - 1
            end = end -1
            gene_seq = seq0.seq[start:end]
        else:
            gene_seq = seq0.seq[start:end].reverse_complement()
        OG_seq[geneid] = gene_seq
    #Write the gene sequences to a fasta file
    output_dir = Path(output_folder)
    file_name = str(output_dir / OG_name) + '.fna'
    with open(file_name, 'w') as f:
        for geneid, gene_seq in OG_seq.items():
            f.write('>' + geneid + '\n' + str(gene_seq) + '\n')
    command = ['clustalw2', '-align', f'-INFILE={file_name}']
    subprocess.run(command, capture_output=True, text=True)
    #Now read the dnd file
    dnd_file = str(output_dir / OG_name) + '.dnd'
    tree = Tree(dnd_file)
    leaves = tree.get_leaf_names()
    distance = set(tree.get_distance(leaf) for leaf in leaves)
    fna0_list = fna_list
    if len(distance) >= 4:
        gff_df_gene['start'] = gff_df_gene['start'] - 200
        gff_df_gene['end'] = gff_df_gene['end'] + 200
        #Generate the new gene sequences
        OG_seq = {}
        for i in range(len(gff_df_gene)):
            line = gff_df_gene.iloc[i]
            seqid = line['seqid']
            start = line['start']
            end = line['end']
            strand = line['strand']
            geneid = line['attributes']            
            seq0 = [x for x in fna0_list if x.id == seqid][0]
            if strand == '+':
                start = start - 1
                end = end -1
                gene_seq = seq0.seq[start:end]
            else:
                gene_seq = seq0.seq[start:end].reverse_complement()
            OG_seq[geneid] = gene_seq

        #Write the gene sequences to a fasta file
        result_dir = Path(result_folder)
        file_name = str(output_dir / OG_name) + '.target.fna'
        output_file = str(result_dir / OG_name) + '.clustalo.aln'
        with open(file_name, 'w') as f:
            for geneid, gene_seq in OG_seq.items():
                f.write('>' + geneid + '\n' + str(gene_seq) + '\n')
        command = ['clustalo', '-i', f'{file_name}', '-o', f'{output_file}']
        subprocess.run(command, capture_output=True, text=True)
        print(f"{OG_name} is candidate")            
    return

def process_files_parallel(OG_list, output_folder, num_threads):
    # Create a pool of worker processes
    pool = mp.Pool(processes=num_threads)
    # Create a partial function with the output_file argument
    func = partial(find_candidate, output_folder=output_folder)
    # Map the function to all files in parallel
    pool.map(func, OG_list)
    # Close the pool and wait for all processes to finish
    pool.close()
    pool.join()

if __name__ == "__main__":
    #Get the path of $PBS_JOBFS
    workdir = os.getenv('PBS_JOBFS')
    workdir = os.path.join(workdir, 'SOG')
    os.chdir(workdir)
    #assign each OG as each task
    OG_list = pd.read_csv("OG_output/Orthogroups_SingleCopyOrthologues.txt", sep='\t', header=None)
    OG_list = OG_list[0].tolist()
    num_threads = 32
    output_folder = os.path.join(workdir, 'tmp')
    result_folder = os.path.join(workdir, 'MSA')
    fasta_folder = os.path.join(workdir, 'genome')
    fasta_files = glob.glob(fasta_folder + '/*.scaffolds.fa')
     #read fasta files and add to seqio_list
    fna_list = []
    for fasta_file in fasta_files:
        fna = list(SeqIO.parse(fasta_file, 'fasta'))
        fna_list.extend(fna)
    gff_folder = os.path.join(workdir, 'gff')
    gff_files = glob.glob(gff_folder + '/*.gff3')
    gff_df = pd.DataFrame()
    for gff_file in gff_files:
        gff = pd.read_csv(gff_file, sep='\t', comment='#', header=None)
        gff.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
        gff_df = pd.concat([gff_df, gff])
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    if not os.path.exists(result_folder):
        os.makedirs(result_folder)
    process_files_parallel(OG_list, output_folder, num_threads)

