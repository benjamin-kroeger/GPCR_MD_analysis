import os
import pathlib
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

def compare_clustering(clustering1_path:str,clustering2_path:str):

    clustering1 = pd.read_csv(clustering1_path)
    clustering1_name = os.path.basename(clustering1_path).split('.')[0]
    clustering1.columns = clustering1.columns[:-1].tolist() + [clustering1_name]

    clustering2 = pd.read_csv(clustering2_path)
    clustering2_name = os.path.basename(clustering2_path).split('.')[0]
    clustering2.columns = clustering2.columns[:-1].tolist() +[clustering2_name]


    merged_df = pd.merge(clustering1,clustering2,on=['0'])

    cross_tab = pd.crosstab(merged_df[clustering1_name], merged_df[clustering2_name])
    sns.heatmap(cross_tab, annot=True, fmt='', cmap='YlGnBu')
    plt.show()
    print(merged_df)


if __name__ == '__main__':

    compare_clustering('output_dir/assinged_clusters_schrö.csv', 'output_dir/assinged_clusters_tmali_max.csv')