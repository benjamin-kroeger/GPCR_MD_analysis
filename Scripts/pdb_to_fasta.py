import os

import pandas as pd
from Bio import SeqIO
import re
def pdb_to_fasta():
    # Replace 'your_input.pdb' with your PDB file and 'pdb_sequences.fasta' with the desired FASTA file name.
    output_fasta = "pdb_sequences.fasta"
    input_pdb_folder = "/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/All_pdb_files"
    records = []

    read_rows = []
    name_pattern = re.compile(r'[A-Z0-9]{3,5}')
    for file in os.listdir(input_pdb_folder):
        name = re.findall(name_pattern,file)[0]
        if name in read_rows:
            continue
        with open(os.path.join(input_pdb_folder, file), 'r') as pdbfile:
            for record in SeqIO.parse(pdbfile, "pdb-seqres"):
                record.id = name
                record.description = ''
                records.append(record)
        read_rows.append(name)

    SeqIO.write(records, output_fasta, "fasta")

def strip_phase_clusters():

    phase_df = pd.read_csv('output_dir/assinged_clusters_schrö.csv')
    phase_df['0'] = phase_df['0'].apply(lambda x: x.split('_')[0])
    phase_df = phase_df.drop_duplicates()
    phase_df.to_csv('output_dir/assinged_clusters_schrö_stripped.csv',index=False)


if __name__ == '__main__':
    pdb_to_fasta()