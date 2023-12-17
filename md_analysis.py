import os
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.manifold import TSNE

from SimilarityAnalyser import SimilarityAnalyser

SCHRODINGER = os.environ.get('SCHRODINGER', '/opt/schrodinger2023-2')


def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--method',type=str,required=True)
    parser.add_argument('--md_frames', type=str, help='Path to the md frames. Exclusive with pdb files')
    parser.add_argument('--pdb_files', type=str, help='Path to a folder with all pdb files. Exclusive with md frames')
    parser.add_argument('--working-dir', type=str, help='Path to a dir where the fiels shall be stored')
    parser.add_argument('--skip_vol_calc', action='store_true', default=False)

    args = parser.parse_args()

    return args


def create_clustering(overlap_matrix_name: str):
    # read and prepare overlap matrix
    similarity_matrix = pd.read_csv(overlap_matrix_name, index_col=0)
    # Append .1 suffix to duplicates just like pandas does implicitly to cols
    new_index = similarity_matrix.index.to_series().reset_index(drop=True)
    duplicate_mask = similarity_matrix.index.duplicated(keep='first')
    new_index[duplicate_mask] = new_index[duplicate_mask] + '.1'
    similarity_matrix.index = new_index

    distance_matrix_values = 1-similarity_matrix.to_numpy()

    sns.clustermap(data=distance_matrix_values, method='average')
    plt.show()

    # make an initial heat mapplot
    ax = sns.heatmap(similarity_matrix.to_numpy(), xticklabels=False, yticklabels=False)
    ax.set(title='Volume overlap heatmap before clustering')
    plt.show()
    # create clusters using hirachial linkage
    linkage_data = linkage(distance_matrix_values, method='average', optimal_ordering=True)
    dendrogram(linkage_data)
    plt.show()
    # find the clusters
    cluster_df = fcluster(linkage_data, t=7, criterion='maxclust')
    assigned_cluster_df = pd.DataFrame(data=cluster_df, columns=['cluster'], index=similarity_matrix.index)
    assigned_cluster_df.sort_values(by=['cluster'], inplace=True)
    # reorder the dataframe
    #rows
    similarity_matrix = similarity_matrix.reindex(index=assigned_cluster_df.index)
    #cols
    similarity_matrix = similarity_matrix[assigned_cluster_df.index]
    # create a new heatmap
    manipulated_similarity_matrix = manipulate_heatmap_data(similarity_matrix.to_numpy(), assigned_cluster_df)
    ax = sns.heatmap(manipulated_similarity_matrix, xticklabels=False, yticklabels=False, square=False)
    ax.set(title='Volume overlap heatmap after clustering')
    plt.show()

    return similarity_matrix, assigned_cluster_df


def manipulate_heatmap_data(similarity_matrix: np.ndarray, clusters: pd.DataFrame):
    cluster_data = clusters['cluster'].reset_index(drop=True)
    change_indices = cluster_data[cluster_data != cluster_data.shift(1)].index.tolist()

    masking_data_row = [[1] * similarity_matrix.shape[0]] * (len(change_indices))
    similarity_matrix = np.insert(similarity_matrix, change_indices, masking_data_row, axis=0)

    masking_data_col = [[1] * similarity_matrix.shape[0]] * (len(change_indices))
    similarity_matrix = np.insert(similarity_matrix, change_indices, np.array(masking_data_col).T, axis=1)

    return similarity_matrix


def create_tsne(ordered_overlap_matrix: pd.DataFrame, cluster_assignments: pd.DataFrame):
    tsne = TSNE(n_components=2, random_state=42)
    tsne_values = tsne.fit_transform(ordered_overlap_matrix)

    plt.figure(figsize=(10, 8))
    for cluster_id in cluster_assignments['cluster'].unique():
        plt.scatter(tsne_values[cluster_assignments['cluster'] == cluster_id, 0],
                    tsne_values[cluster_assignments['cluster'] == cluster_id, 1],
                    label=f'cluster_{cluster_id}'
                    )
    plt.title('t-SNE Plot on similarity matrix with Cluster Colors')
    plt.legend()
    plt.show()


def main(args):
    # just one of the ways can be used to generate a similarity matrix
    plt.rcParams["figure.figsize"] = (16, 12)
    extractor = SimilarityAnalyser()

    if args.method == 'schrödinger':
        assert args.md_frames is not None, "Md frames is missing"
        volume_overlap_filename = extractor.compute_similarity_volume_overlap(working_dir=args.working_dir, md_frames=args.md_frames)
        suffix = 'schrö'
    if args.method == 'tm_align':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similarity_tmalign(pdb_files=args.pdb_files)
        suffix = 'tmali'
    if args.method == 'gogo':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similarity_go(pdb_file_dir=args.pdb_files)
        suffix = 'GOGO'
    if args.method == 'vmd':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similaritiy_vmd(pdb_files=args.pdb_files)
        suffix = 'vmd'


    print('Start plotting and clustering')
    similarity_matrix, assigned_clusters = create_clustering(overlap_matrix_name=volume_overlap_filename)
    assigned_clusters.to_csv(f'output_dir/assinged_clusters_{suffix}.csv')
    create_tsne(ordered_overlap_matrix=similarity_matrix, cluster_assignments=assigned_clusters)


if __name__ == '__main__':
    args = parse_args()
    main(args)
