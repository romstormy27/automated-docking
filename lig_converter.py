# import necessary library
import subprocess
import os
from tqdm import tqdm

# dir to mol2 files
ligand_dir = "dataset/processed/mol2"

# loop through all mol2 files and convert into pdbqt for docking readiness
for name in tqdm(os.listdir(ligand_dir)):
    
    output_name = name.replace(".mol2", ".pdbqt")

    subprocess.run(["bash", "ADFRsuite-1.0/bin/run.sh", name, output_name])