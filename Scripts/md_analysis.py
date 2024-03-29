import os
from argparse import ArgumentParser

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from SimilarityAnalyser import SimilarityAnalyser
import fastcluster
SCHRODINGER = os.environ.get('SCHRODINGER', '/opt/schrodinger2023-2')


def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--method', type=str, required=True)
    parser.add_argument('--md_frames', type=str, help='Path to the md frames. Exclusive with pdb files')
    parser.add_argument('--pdb_files', type=str, help='Path to a folder with all pdb files. Exclusive with md frames')
    parser.add_argument('--working-dir', type=str, help='Path to a dir where the fiels shall be stored', required=False)
    parser.add_argument('--tm_align_path', type=str, required=False,
                        default="/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/TMalign")
    parser.add_argument('--output_dir', required=False, type=str, default="output_dir")
    parser.add_argument('--gogo_path', type=str, required=False, default="/home/benjaminkroeger/Downloads/GOGO_master")
    parser.add_argument("--num_clusters", type=int, default=7,
                        help="The number of clusters to use for clustering (default:7). Does not impact similarity computation")
    args = parser.parse_args()

    return args

def assign_clusters_linkage(distance_matrix_values, org_index)->pd.DataFrame:
    # create clusters using hirachial linkage
    # linkage_data = linkage(distance_matrix_values, method='average', optimal_ordering=True)
    linkage_data = fastcluster.linkage(distance_matrix_values, method='average', metric='euclidean')
    dendrogram(linkage_data)
    plt.show()
    # find the clusters
    # cluster_df = fcluster(linkage_data, t=num_clusters, criterion='maxclust')
    cluster_df = fcluster(linkage_data, t=3, criterion='distance')
    assigned_cluster_df = pd.DataFrame(data=cluster_df, columns=['cluster'], index=org_index)
    assigned_cluster_df.sort_values(by=['cluster'], inplace=True)

    return assigned_cluster_df


def assign_clusters_kmeans(distance_matrix_values, org_index, num_clusters=7):
    # Perform K-means clustering
    kmeans_model = KMeans(n_clusters=num_clusters, random_state=42)
    clusters = kmeans_model.fit_predict(distance_matrix_values)

    assigned_cluster_df = pd.DataFrame(data=clusters + 1, columns=['cluster'], index=org_index)
    assigned_cluster_df.sort_values(by=['cluster'], inplace=True)

    return assigned_cluster_df

def create_clustering(overlap_matrix_name: str,num_clusters:int):
    # read and prepare overlap matrix
    similarity_matrix = pd.read_csv(overlap_matrix_name, index_col=0)
    # Append .1 suffix to duplicates just like pandas does implicitly to cols
    new_index = similarity_matrix.index.to_series().reset_index(drop=True)
    duplicate_mask = similarity_matrix.index.duplicated(keep='first')
    new_index[duplicate_mask] = new_index[duplicate_mask] + '.1'
    similarity_matrix.index = new_index

    distance_matrix_values = 1 - similarity_matrix.to_numpy()

    sns.clustermap(data=distance_matrix_values, method='average')
    plt.show()

    # make an initial heat mapplot
    ax = sns.heatmap(similarity_matrix.to_numpy(), xticklabels=False, yticklabels=False)
    ax.set(title='Volume overlap heatmap before clustering')
    plt.show()
    #assigned_cluster_df = assign_clusters_linkage(distance_matrix_values=distance_matrix_values, org_index=similarity_matrix.index)
    assigned_cluster_df = assign_clusters_kmeans(distance_matrix_values=distance_matrix_values, org_index=similarity_matrix.index)
    # reorder the dataframe
    # rows
    similarity_matrix = similarity_matrix.reindex(index=assigned_cluster_df.index)
    # cols
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
    extractor = SimilarityAnalyser(tm_align_path=args.tm_align_path,
                                   gogo_dir_path=args.gogo_path,
                                   output_dir=args.output_dir)

    if args.method == 'schrödinger':
        assert args.md_frames is not None, "Md frames is missing"
        volume_overlap_filename = extractor.compute_similarity_volume_overlap(working_dir=args.working_dir, md_frames=args.md_frames)
        suffix = 'schrö'
    if args.method == 'tm_align':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similarity_tmalign(pdb_files=args.pdb_files)
        suffix = 'tm_align'
    if args.method == 'gogo':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similarity_go(pdb_files=args.pdb_files)
        suffix = 'GOGO'
    if args.method == 'vmd':
        assert args.pdb_files is not None, "Pdb file path is missing"
        volume_overlap_filename = extractor.compute_similaritiy_vmd(pdb_files=args.pdb_files)
        suffix = 'vmd'

    print('Start plotting and clustering')
    similarity_matrix, assigned_clusters = create_clustering(overlap_matrix_name=volume_overlap_filename,num_clusters=args.num_clusters)
    assigned_clusters.to_csv(f'{args.output_dir}/assinged_clusters_{suffix}.csv')
    create_tsne(ordered_overlap_matrix=similarity_matrix, cluster_assignments=assigned_clusters)


if __name__ == '__main__':
    args = parse_args()
    main(args)
