import copy
import csv
import os
import re
import shutil
import subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime

import pandas as pd
from tqdm import tqdm

SCHRODINGER = os.environ.get('SCHRODINGER', '/opt/schrodinger2023-2')


class SimilarityAnalyser:

    def __init__(self):
        self.tmalign_path = "/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/TMalign"

    def compute_similarity_tmalign(self, pdb_files: str) -> str:
        triplet_output_path = 'output_dir/output_triplets_max.csv'
        similarity_output_matrix_path = 'output_dir/output_tmalign_simliarity_max.csv'
        if os.path.exists(similarity_output_matrix_path):
            return similarity_output_matrix_path

        assert os.path.isdir(pdb_files), "The directory to the pdb files does not exist"
        files = os.listdir(pdb_files)
        file_sparse = copy.deepcopy(files)
        result_list = []
        jobs = []
        if os.path.exists(triplet_output_path):
            triplet_df = pd.read_csv(triplet_output_path, header=None)
        else:
            with ProcessPoolExecutor(max_workers=16) as executor:
                for pdb_file1 in files:

                    for pdb_file2 in file_sparse:
                        job = executor.submit(self._get_tm_align_score, os.path.join(pdb_files, pdb_file1), os.path.join(pdb_files, pdb_file2))
                        jobs.append(job)

                    file_sparse.remove(pdb_file1)

            print('All jobs submitted')

            with tqdm(total=len(files) // 2) as pbar:
                with open(triplet_output_path, 'w') as f:
                    writer = csv.writer(f)
                    for finished_job in as_completed(jobs):
                        result = finished_job.result()
                        result_list.append(result)
                        writer.writerow(result)
                        pbar.update(1)

            triplet_df = pd.DataFrame(data=result_list, index=None, columns=None)

        similarity_matrix = self._convert_tripelts_into_matrix(triplet_df)
        similarity_matrix.to_csv(similarity_output_matrix_path)
        return similarity_output_matrix_path

    def _convert_tripelts_into_matrix(self, triplet_df: pd.DataFrame) -> pd.DataFrame:
        assert triplet_df[0].nunique() == triplet_df[1].nunique()
        colnames = triplet_df[0].unique()
        colnames.sort()

        triplet_df = triplet_df.pivot(index=0, columns=1, values=2)
        triplet_df = triplet_df.reindex(columns=colnames, index=colnames)

        similarity_matrix = triplet_df.add(triplet_df[~triplet_df.isna()].T, fill_value=0)
        return similarity_matrix

    def _get_tm_align_score(self, path_to_pdb1, path_to_pdb2) -> tuple[str, str, float]:
        tm_align_out = (self._run_cmdline(' '.join([self.tmalign_path, path_to_pdb1, path_to_pdb2])))
        tm_align_out = tm_align_out.decode('utf-8')

        return os.path.basename(path_to_pdb1).rstrip('.pdb'), os.path.basename(path_to_pdb2).rstrip('.pdb'), self._process_tm_align_output(tm_align_out)

    def _process_tm_align_output(self, tm_align_output: str) -> float:
        score_pattern = re.compile(r'TM-score= (\d+\.\d+)')
        scores = [float(x.group(1)) for x in re.finditer(score_pattern, tm_align_output)]

        return max(scores)

    def compute_similarity_volume_overlap(self, working_dir: str, md_frames: str) -> str:
        similarity_matrix_path = 'output_dir/output_phase_vol.csv'
        if os.path.exists(similarity_matrix_path):
            return similarity_matrix_path

        # create new dir if it does not exist
        if not os.path.isdir(working_dir):
            os.mkdir(working_dir)
        else:
            os.mkdir(working_dir + '_' + datetime.now().strftime('%Y-%m-%d-%H:%M'))

        # copy the file if not present
        expected_file_loc = os.path.join(working_dir, os.path.basename(md_frames))
        if not os.path.isfile(expected_file_loc):
            shutil.copy(md_frames, expected_file_loc)

        # run the analysis
        self._run_volume_cluster(os.path.basename(md_frames), working_dir)
        volume_overlap_filename = self._run_phase_volCalc(os.path.basename(md_frames), working_dir)

        similarity_matrix = pd.read_csv(volume_overlap_filename)
        similarity_matrix.index = similarity_matrix['ID']
        similarity_matrix.drop(columns=['ID'], inplace=True)
        similarity_matrix.index.name = 0
        similarity_matrix.to_csv(similarity_matrix_path)

        return similarity_matrix_path

    def _run_cmdline(self, cmd: str) -> bytes:
        process = subprocess.Popen(cmd,
                                   stderr=subprocess.PIPE,
                                   stdout=subprocess.PIPE,
                                   stdin=subprocess.PIPE,
                                   shell=True
                                   )

        stdout, stderr = process.communicate(timeout=None)

        # print(stdout.decode('utf-8'))
        # print(stderr.decode('utf-8'))

        return stdout

    def _run_volume_cluster(self, file_name: str, path_to_working_dir):
        volume_cluster_cmd = f'cd {path_to_working_dir}; ' \
                             f'{SCHRODINGER}{os.sep}run volume_cluster.py -a backbone -g 1 -j {file_name.split(".")[0]} -l Average -sc -WAIT -HOST localhost {file_name}'
        print('Starting volume_clust')
        self._run_cmdline(volume_cluster_cmd)

    def _run_phase_volCalc(self, file_name: str, path_to_working_dir: str) -> str:
        phase_volCalc_cmd = f'cd {path_to_working_dir}; ' \
                            f'{SCHRODINGER}/utilities/phase_volCalc -mae1 asl.mae -mae2 asl.mae -out {file_name.split(".")[0]}_titles.csv -grid 1.0 -titles1 -titles2 -scores'
        print('Starting phase_volCalc')
        self._run_cmdline(phase_volCalc_cmd)

        return f'{file_name.split(".")[0]}_titles.csv'
