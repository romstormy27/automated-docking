import os
import subprocess
from tqdm import tqdm

rec_dir = "dataset/raw/pdb"

def pdb_to_pdbqt(input_path:str, output_path:str):

    for name in tqdm(os.listdir(rec_dir)):
    
        output_name = name.replace(".pdb", ".pdbqt")

        subprocess.run(["./ADFRsuite-1.0/bin/prepare_receptor", "-r", input_path+name, "-o", output_path+output_name])