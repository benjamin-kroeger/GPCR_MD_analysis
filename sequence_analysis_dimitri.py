import os

from Bio import SeqIO, Align
from Bio.Blast import Applications
from Bio.Align.Applications import ClustalOmegaCommandline

os.environ['PATH'] = f"/home/benjaminkroeger/Downloads/clustalo:{os.environ['PATH']}"

def create_msa(sequences):


    # Create a MuscleCommandline object to run the MUSCLE alignment tool.
    clustalw = ClustalOmegaCommandline("/home/benjaminkroeger/Downloads/clustalo",infile="pdb_sequences.fasta", outfile="output_alig.fasta",wrap=300)
    clustalw()

    # You can also use other alignment methods like Clustal Omega or T-Coffee with appropriate command lines.

    # Read the MSA result from the output file (replace 'pdb_sequences.fasta' with your output file).
    aligned_sequences = list(SeqIO.parse("pdb_sequences.fasta", "fasta"))

    # Print or use the aligned sequences as needed.
    for record in aligned_sequences:
        print(record.id)
        print(record.seq)


if __name__ == "__main__":
    sequences = list(SeqIO.parse("pdb_sequences.fasta", "fasta"))
    create_msa(sequences)