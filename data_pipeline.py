# necessary for compound sdf acquisition
import pubchempy as pcp
from openbabel import openbabel
import pandas as pd

# necessary for preparing ligand and receptor
import subprocess
import os
from tqdm import tqdm
import os

# general import
import yaml
from utils import load_config

# ===============================LIGAND PREPARATION============================

def get_sdf_from_pubchem(cid_list:list, raw_sdf_path:str):

    '''
    all you have to do is prepare the list of ligand identifier as a python list
    and loop through them to get the sdf files
    
    '''
    print("Downloading...")
    for ligand in tqdm(set(cid_list)):

        try:
            # download sdf from pubchem
            pcp.download('SDF', raw_sdf_path+str(ligand)+".sdf", ligand, "cid", overwrite=True)

        except pcp.NotFoundError:

            print(ligand, "is not found")
    print()

def convert_sdf_to_mol2(sdf_path:str, mol2_path:str):

    # loop through all sdf files to convert into mol2
    print("Converting(sdf->mol2)...")
    for sdf in tqdm(os.listdir(sdf_path)):

        # file instance
        ob = openbabel.OBConversion()

        # set in out formats
        ob.SetInAndOutFormats("sdf", "mol2")

        # instantiate openbabel object
        mol = openbabel.OBMol()

        # read sdf file and add to openbabel object
        ob.ReadFile(mol, sdf_path+sdf)

        output_name = sdf.replace(".sdf", ".mol2")

        # write open babel object as mol2 file
        ob.WriteFile(mol, mol2_path+output_name)
    print()

def convert_mol2_to_pdbqt(mol2_path:str):

    # loop through all mol2 files and convert into pdbqt for docking readiness
    print("Converting(mol2->pdbqt)...")
    for name in tqdm(os.listdir(mol2_path)):
        
        # define output format
        output_name = name.replace(".mol2", ".pdbqt")

        # run shell script
        subprocess.run(["bash", "ADFRsuite-1.0/bin/run.sh", name, output_name])
    print()

# ===============================PROTEIN PREPARATION============================

def convert_pdb_to_pdbqt(pdb_path:str, pdbqt_path:str):

    print("Converting(pdb->pdbqt)...")
    for name in tqdm(os.listdir(pdb_path)):
    
        output_name = name.replace(".pdb", ".pdbqt")

        subprocess.run(["./ADFRsuite-1.0/bin/prepare_receptor", "-r", 
                        pdb_path+name, "-o", pdbqt_path+output_name])
    print()

# ===============================DRIVER CODE============================

if __name__=="__main__":

    config = load_config()

    cid_list = config["ligand_cid_list"]

    # downloaded sdf path
    raw_sdf_path = config["raw_sdf_path"]
    # mol2 path
    mol2_path = config["mol2_path"]
    
    # receptor pdb path
    receptor_pdb_path = config["receptor_pdb_path"]
    # receptor pdbtq path
    receptor_pdbqt_path = config["receptor_pdbqt_path"]

    # run get sdf from pubchem function
    get_sdf_from_pubchem(cid_list=cid_list, raw_sdf_path=raw_sdf_path)

    # run convert sdf to mol2 function
    convert_sdf_to_mol2(sdf_path=raw_sdf_path, mol2_path=mol2_path)

    # run convert mol2 to sdf function
    convert_mol2_to_pdbqt(mol2_path=mol2_path)

    # run convert recepter pdb to pdbqt
    convert_pdb_to_pdbqt(pdb_path=receptor_pdb_path, pdbqt_path=receptor_pdbqt_path)

