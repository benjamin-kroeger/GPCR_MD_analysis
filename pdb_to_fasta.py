import os

from Bio import SeqIO

# Replace 'your_input.pdb' with your PDB file and 'pdb_sequences.fasta' with the desired FASTA file name.
output_fasta = "pdb_sequences.fasta"
input_pdb_folder = "/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/All_pdb_files"
records = []

read_rows = []
for file in os.listdir(input_pdb_folder):
    name = file.split('_')[0]
    if name in read_rows:
        continue
    with open(os.path.join(input_pdb_folder, file), 'r') as pdbfile:
        for record in SeqIO.parse(pdbfile, "pdb-seqres"):
            record.id = file.rstrip('.pdb').split('_')[0]
            record.description = ''
            records.append(record)
    read_rows.append(name)

SeqIO.write(records, output_fasta, "fasta")
