import csv
import os

import MDAnalysis.analysis.rms
import dendropy
import numpy as np
from Bio import SeqIO, pairwise2
from Bio.Align.Applications import ClustalOmegaCommandline
from sklearn.cluster import KMeans

os.environ['PATH'] = f"/home/benjaminkroeger/Downloads/clustalo:{os.environ['PATH']}"


def calculate_rmsd_sim_pairs(file_path1, file_path2):
    file1 = MDAnalysis.Universe(file_path1)
    file2 = MDAnalysis.Universe(file_path2)
    R = MDAnalysis.analysis.rms.RMSD(file1, file2, select="backbone")
    R.run()
    rmsd = R.rmsd[2]

    sequences = []
    for pdb_file in [file_path1, file_path2]:
        with open(pdb_file, "r") as handle:
            for record in SeqIO.parse(handle, "pdb-atom"):
                sequences.append(str(record.seq))

        # Perform pairwise sequence alignment
    alignments = pairwise2.align.globalxx(sequences[0], sequences[1])

    # Calculate sequence identity
    alignment = alignments[0]
    seq_len = max(len(sequences[0]), len(sequences[1]))
    identity = alignment[2] / seq_len * 100  # Calculating percentage identity

    return (rmsd,identity)




def create_msa():
    # Create a MuscleCommandline object to run the MUSCLE alignment tool.
    clustalw = ClustalOmegaCommandline("/home/benjaminkroeger/Downloads/clustalo", infile="Alignments_and_fasta/pdb_sequences.fasta",
                                       outfile="Alignments_and_fasta/clustalO_MSA.fasta", wrap=300,
                                       clusteringout='cluster.txt', force=True,
                                       outfmt="clu", guidetree_out="Alignments_and_fasta/tree.out", outputorder="tree-order",
                                       seqtype='protein')
    clustalw()


def apply_dendropy(num_clusters):
    with open('Alignments_and_fasta/tree.out', 'r') as f:
        tree_str = ''.join(f.readlines())

    # Read the Newick tree data
    tree = dendropy.Tree.get(data=tree_str, schema="newick")

    # Extract the leaf labels (sequence names) and compute a feature matrix
    leaf_labels = [leaf.taxon.label for leaf in tree.leaf_node_iter()]
    feature_matrix = np.zeros((len(leaf_labels), 1))  # 1 feature dimension

    # Perform k-means clustering
    kmeans = KMeans(n_clusters=num_clusters, random_state=0).fit(feature_matrix)

    # Assign each sequence to a cluster
    clusters = {}
    for i, label in enumerate(leaf_labels):
        cluster_label = kmeans.labels_[i]
        if cluster_label not in clusters:
            clusters[cluster_label] = []
        clusters[cluster_label].append(label)

    # Save the cluster assignments to a CSV file
    with open('clustalo_clusters.csv', 'w', newline='') as f:
        csv_writer = csv.writer(f)
        csv_writer.writerow(['cluster', 'sequence'])
        for cluster_label, cluster_members in clusters.items():
            for seq in cluster_members:
                csv_writer.writerow([cluster_label, seq])


if __name__ == "__main__":
    create_msa()
