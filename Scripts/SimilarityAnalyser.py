import copy
import csv
import os
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed, wait
from datetime import datetime
from io import StringIO

import numpy as np
import requests
import pandas as pd
from Bio import SeqIO
from tqdm import tqdm

SCHRODINGER = os.environ.get('SCHRODINGER', '/opt/schrodinger2023-2')
DENSITY_MAP_TCL = """
mol new {0} type pdb
set atomselect{1} [atomselect {1} "all"]
volmap density atomselect{1} -res 1.0 -radscale 1.0 -weight mass -o {2}
"""
VOLTOOL_TCL ="""
voltool correlate -i1 {0} -i2 {1}
"""

class SimilarityAnalyser:

    def __init__(self):
        self.tmalign_path = "/home/benjaminkroeger/Documents/Master/Master_2_Semester/Internship/TMalign"
        self.gogo_path = "/home/benjaminkroeger/Downloads/GOGO_master"
        self.output_dir = "output_dir"

        if not os.path.isdir(self.output_dir):
            os.mkdir(self.output_dir)

    def compute_similaritiy_vmd(self, pdb_files: str) -> str:
        triplet_output_path = os.path.join(self.output_dir,'output_triplets_vmd.csv')
        similarity_output_matrix_path = os.path.join(self.output_dir,'output_vmd_simliarity.csv')

        assert os.path.isdir(pdb_files), "The directory to the pdb files does not exist"
        tmp_dir = tempfile.TemporaryDirectory()
        print(tmp_dir.name)
        files = os.listdir(pdb_files)

        if os.path.exists(triplet_output_path):
            triplet_df = pd.read_csv(triplet_output_path, header=None)
        else:

            whole_path_pdb_files = [os.path.join(pdb_files, x) for x in files]
            density_jobs = []
            with ProcessPoolExecutor(max_workers=16) as executor:
                for i in range(0, len(whole_path_pdb_files), 100):
                    density_jobs.append(executor.submit(self._create_density_maps, whole_path_pdb_files[i:i + 100], tmp_dir))

            wait(density_jobs)
            print('Created density maps for each pdb file')
            pairings = []
            pdb_ids = [x.rstrip('.pdb') for x in files]
            pdb_ids_sparse = copy.deepcopy(pdb_ids)
            for id1 in pdb_ids:
                for id2 in pdb_ids_sparse:
                    pairings.append((id1,id2))

                pdb_ids_sparse.remove(id1)

            voltool_jobs = []
            with ProcessPoolExecutor(max_workers=16) as executor:
                for i in range(0,len(pairings),400):

                    job = executor.submit(self._compute_ccc_from_density,pairings[i:i+400],tmp_dir)
                    voltool_jobs.append(job)

            result_list = []
            with open(triplet_output_path, 'w') as f:
                writer = csv.writer(f)
                for finished_job in as_completed(voltool_jobs):
                    result = finished_job.result()
                    result_list.extend(result)
                    for entry in result:
                        writer.writerow(entry)

            triplet_df = pd.DataFrame(data=result_list, index=None, columns=None)

        tmp_dir.cleanup()

        similarity_matrix = self._convert_tripelts_into_matrix(triplet_df)
        similarity_matrix.to_csv(similarity_output_matrix_path)

        return similarity_output_matrix_path

    def _create_density_maps(self, pdb_files: list[str], tmp_dir):
        # create density script location and file name
        density_tcl_script = ''
        density_tcl_script_hash = abs(hash(''.join(pdb_files)))
        density_tcl_script_path = os.path.join(tmp_dir.name, f'density_script_{density_tcl_script_hash}.tcl')

        # create the script
        for i, pdb_file_path in enumerate(pdb_files):
            density_file_path = os.path.join(tmp_dir.name, os.path.basename(pdb_file_path).rstrip('pdb') + 'dx')
            density_tcl_script = density_tcl_script + DENSITY_MAP_TCL.format(pdb_file_path, i, density_file_path)
        # write the script
        with open(density_tcl_script_path, 'w') as f:
            f.write(density_tcl_script)

        # execute the script
        output_path = os.path.join(tmp_dir.name, f"dummy_{density_tcl_script_hash}.log")
        vmd_cmd = f"vmd -dispdev text -eofexit < {density_tcl_script_path} > {output_path}"
        self._run_cmdline(vmd_cmd)

    def _compute_ccc_from_density(self,pairings:list[tuple[str,str]], tmpdir) -> list[tuple[str,str,float]]:

        voltool_tcl = ''
        voltool_script_hash = abs(hash(str(pairings)))
        voltool_script_path = os.path.join(tmpdir.name,f'voltool_script_{voltool_script_hash}')
        voltool_script_output_path = os.path.join(tmpdir.name,f'voltool_output_{voltool_script_hash}')
        for id1,id2 in pairings:
            path_density_1 = os.path.join(tmpdir.name, f'{id1}.dx')
            path_density_2 = os.path.join(tmpdir.name, f'{id2}.dx')

            voltool_tcl = voltool_tcl + VOLTOOL_TCL.format(path_density_1,path_density_2)

        with open(voltool_script_path,'w') as f:
            f.write(voltool_tcl)

        vmd_cmd = f"vmd -dispdev text -eofexit < {voltool_script_path} > {voltool_script_output_path}"
        self._run_cmdline(vmd_cmd)
        scores = self._parse_vmd_output(voltool_script_output_path)

        return [(pair[0],pair[1],score) for pair,score in zip(pairings,scores)]

    def _parse_vmd_output(self, vmd_output_path: str) -> list[float]:
        with open(vmd_output_path, 'r') as vmd_output_handle:
            output = ''.join(vmd_output_handle.readlines())
            scores = re.findall(r'^\d\.\d+',output,re.MULTILINE)

            return [float(x) for x in scores]

    def compute_similarity_tmalign(self, pdb_files: str) -> str:
        triplet_output_path = os.path.join(self.output_dir,'output_triplets_min.csv')
        similarity_output_matrix_path = os.path.join(self.output_dir,'output_tmalign_simliarity_min.csv')
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
        assert triplet_df[triplet_df.columns[0]].nunique() == triplet_df[triplet_df.columns[1]].nunique()
        colnames = triplet_df[triplet_df.columns[0]].unique()
        colnames.sort()

        triplet_df = triplet_df.pivot(index=triplet_df.columns[0], columns=triplet_df.columns[1], values=triplet_df.columns[2])
        triplet_df = triplet_df.reindex(columns=colnames, index=colnames)

        similarity_matrix = triplet_df.add(triplet_df[~triplet_df.isna()].T, fill_value=0)

        if similarity_matrix.iloc[0,0] > 1:
            for i in range(min(similarity_matrix.shape)):
                similarity_matrix.iat[i, i] = 1

        return similarity_matrix

    def compute_similarity_go(self, pdb_files: str):
        similarity_output_matrix_path = os.path.join(self.output_dir,'output_gogo_simliarity.csv')
        pdb_ids = set()
        name_pattern = re.compile(r'[A-Z0-9]{3,5}')
        for file in os.listdir(pdb_files):
            name = re.findall(name_pattern, file)[0]
            pdb_ids.add(name)

        print('Getting go annotations')
        go_dict = {pdb_id: self._get_go_annotations(self._convert_pdb_to_uniprot(pdb_id)) for pdb_id in pdb_ids}
        go_dict = dict(sorted(go_dict.items()))
        print('writing inputs')
        gogo_input_file = 'gogo_input.txt'
        with open(gogo_input_file, 'w') as gogo_input:
            for pdb_id, gos in go_dict.items():
                gogo_input.write(f'{pdb_id} {" ".join(gos)}\n')

        output_file_path = self._run_gogo(gogo_input_file)
        gogo_output_triplets = self._parse_gogo_output(output_file_path)

        gogo_sim_matrix = self._convert_tripelts_into_matrix(gogo_output_triplets)
        for i in range(min(gogo_sim_matrix.shape)):
            gogo_sim_matrix.iat[i, i] = 1
        gogo_sim_matrix.dropna(axis=0, thresh=len(gogo_sim_matrix) // 2, inplace=True)
        gogo_sim_matrix.dropna(axis=1, thresh=len(gogo_sim_matrix) // 2, inplace=True)
        gogo_sim_matrix.to_csv(similarity_output_matrix_path)

        return similarity_output_matrix_path

    def _run_gogo(self, input_file):
        cwd = os.getcwd()
        gogo_output = os.path.join(cwd, "gogo_output.tsv")
        gogo_cmd = (f'cd {self.gogo_path}; '
                    f'perl gene_list_comb.pl {os.path.join(cwd, input_file)} {gogo_output} {os.path.join(cwd, "cluster.txt")}')
        print('starting gogo run')
        self._run_cmdline(gogo_cmd)
        return gogo_output

    def _convert_pdb_to_uniprot(self, pdb_id: str) -> str:
        api_url = f"https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%28xref%3A{pdb_id}%29%29"

        resp = requests.get(api_url)
        resp_fasta = resp.text
        uniprot_id = resp_fasta.split('\n')[0].split('|')[1]

        return uniprot_id

    def _get_go_annotations(self, uniprot_id) -> list[str]:

        swiss_entry_url = f"https://rest.uniprot.org/uniprotkb/stream?format=txt&query=%28%28accession%3A{uniprot_id}%29%29"
        swiss_entry = requests.get(swiss_entry_url)

        swiss_rec = next(SeqIO.parse(StringIO(swiss_entry.text), 'swiss'))
        gos = [x.split(':', maxsplit=1)[1] for x in swiss_rec.dbxrefs if x.startswith('GO')]

        return gos

    def _parse_gogo_output(self, output_file_path):
        gogo_output_df = pd.read_csv(output_file_path, sep=' ', names=['GPCR1', 'GPCR2', '_BPO', 'BPO', '_CCO', 'CCO', '_MFO', 'MFO'])
        gogo_output_df.drop(columns=['BPO', 'CCO', '_BPO', '_CCO', '_MFO'], inplace=True)
        return gogo_output_df

    def _get_tm_align_score(self, path_to_pdb1, path_to_pdb2) -> tuple[str, str, float]:
        tm_align_out = (self._run_cmdline(' '.join([self.tmalign_path, path_to_pdb1, path_to_pdb2])))
        tm_align_out = tm_align_out.decode('utf-8')

        return os.path.basename(path_to_pdb1).rstrip('.pdb'), os.path.basename(path_to_pdb2).rstrip('.pdb'), self._process_tm_align_output(
            tm_align_out)

    def _process_tm_align_output(self, tm_align_output: str) -> float:
        score_pattern = re.compile(r'TM-score= (\d+\.\d+)')
        scores = [float(x.group(1)) for x in re.finditer(score_pattern, tm_align_output)]

        return min(scores)

    def compute_similarity_volume_overlap(self, working_dir: str, md_frames: str) -> str:
        similarity_matrix_path = os.path.join(self.output_dir,'output_phase_vol.csv')
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

        print(stdout.decode('utf-8'))
        print(stderr.decode('utf-8'))

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
