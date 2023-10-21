import os
import subprocess
from argparse import ArgumentParser
import shutil
from datetime import datetime
from glob import glob

import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt

import seaborn as sns

import pandas as pd
import sys

SCHRODINGER = os.environ.get('SCHRODINGER', '/opt/schrodinger2023-2')


def parse_args():
    parser = ArgumentParser()

    parser.add_argument('--md-frames', type=str, help='Path to the md frames', required=True)
    parser.add_argument('--working-dir', type=str, help='Path to a dir where the fiels shall be stored', required=True)
    parser.add_argument('--skip_vol_calc', action='store_true',default=False)

    args = parser.parse_args()

    return args


def _run_cmdline(cmd: str):
    process = subprocess.Popen(cmd,
                               stderr=subprocess.PIPE,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               shell=True
                               )

    stdout, stderr = process.communicate(timeout=None)

    print(stdout.decode('utf-8'))
    print(stderr.decode('utf-8'))


def run_volume_cluster(file_name: str, path_to_working_dir):
    volume_cluster_cmd = f'cd {path_to_working_dir}; ' \
                         f'{SCHRODINGER}{os.sep}run volume_cluster.py -a backbone -g 1 -j {file_name.split(".")[0]} -l Average -sc -WAIT -HOST localhost {file_name}'
    print('Starting volume_clust')
    _run_cmdline(volume_cluster_cmd)


def run_phase_volCalc(file_name: str, path_to_working_dir: str)-> str:
    phase_volCalc_cmd = f'cd {path_to_working_dir}; ' \
                        f'{SCHRODINGER}/utilities/phase_volCalc -mae1 asl.mae -mae2 asl.mae -out {file_name.split(".")[0]}_titles.csv -grid 1.0 -titles1 -titles2 -scores'
    print('Starting phase_volCalc')
    _run_cmdline(phase_volCalc_cmd)

    return f'{file_name.split(".")[0]}_titles.csv'



def create_clustering(overlap_matrix_name: str, path_to_working_dir: str):
    # read and prepare overlap matrix
    similarity_matrix = pd.read_csv(os.path.join(path_to_working_dir, overlap_matrix_name))
    similarity_matrix.index = similarity_matrix['ID']
    similarity_matrix.drop(columns=['ID'], inplace=True)
    similarity_matrix_values = similarity_matrix.to_numpy()


    sns.clustermap(data=similarity_matrix_values,method='average')
    plt.show()







    # make an initial heat mapplot
    ax = sns.heatmap(similarity_matrix.to_numpy(),xticklabels=False,yticklabels=False)
    ax.set(title='Volume overlap heatmap before clustering')
    plt.show()
    # create clusters using hirachial linkage
    linkage_data = linkage(similarity_matrix_values, method='average', optimal_ordering=True)
    dendrogram(linkage_data)
    plt.show()
    # find the clusters
    cluster_df = fcluster(linkage_data, t=7, criterion='maxclust')
    assigned_cluster_df = pd.DataFrame(data=cluster_df, columns=['cluster'], index=similarity_matrix.index)
    assigned_cluster_df.sort_values(by=['cluster'], inplace=True)
    # reorder the dataframe
    similarity_matrix = similarity_matrix.reindex(index=assigned_cluster_df.index)
    similarity_matrix = similarity_matrix[assigned_cluster_df.index]
    # create a new heatmap
    manipulated_similarity_matrix = manipulate_heatmap_data(similarity_matrix.to_numpy(),assigned_cluster_df)
    ax = sns.heatmap(manipulated_similarity_matrix,xticklabels=False,yticklabels=False,square=False)
    ax.set(title='Volume overlap heatmap after clustering')
    plt.show()



    return similarity_matrix, assigned_cluster_df

def manipulate_heatmap_data(similarity_matrix:np.ndarray, clusters:pd.DataFrame):

    cluster_data = clusters['cluster'].reset_index(drop=True)
    change_indices = cluster_data[cluster_data != cluster_data.shift(1)].index.tolist()

    masking_data_row = [[1]*similarity_matrix.shape[0]]*(len(change_indices))
    similarity_matrix = np.insert(similarity_matrix,change_indices,masking_data_row,axis=0)

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
    plt.rcParams["figure.figsize"] = (16, 12)

    if not args.skip_vol_calc:
        # create new dir if it does not exist
        if not os.path.isdir(args.working_dir):
            os.mkdir(args.working_dir)
        else:
            os.mkdir(args.working_dir +'_'+ datetime.now().strftime('%Y-%m-%d-%H:%M'))

        # copy the file if not present
        expected_file_loc = os.path.join(args.working_dir, os.path.basename(args.md_frames))
        if not os.path.isfile(expected_file_loc):
            shutil.copy(args.md_frames, expected_file_loc)

        # run the analysis
        run_volume_cluster(os.path.basename(args.md_frames), args.working_dir)
        volume_overlap_filename = run_phase_volCalc(os.path.basename(args.md_frames), args.working_dir)
    else:
        volume_overlap_filename = os.path.basename(glob(args.working_dir + '/*titles.csv')[0])
    print('Start plotting and clustering')
    similarity_matrix, assigned_clusters = create_clustering(overlap_matrix_name=volume_overlap_filename,path_to_working_dir=args.working_dir)
    create_tsne(ordered_overlap_matrix=similarity_matrix, cluster_assignments=assigned_clusters)

if __name__ == '__main__':
    args = parse_args()
    main(args)

