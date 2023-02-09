from openbabel import openbabel
from openbabel import pybel
from tqdm import tqdm
import os

def convert_pdbqt_to_sdf(pdbqt_path:str, sdf_path:str):

    # file instance
    ob = openbabel.OBConversion()

    # set in out formats
    ob.SetInAndOutFormats("pdbqt", "sdf")

    # instantiate openbabel object
    mol = openbabel.OBMol()

    # read sdf file and add to openbabel object
    ob.ReadFile(mol, pdbqt_path)

    # write open babel object as mol2 file
    ob.WriteFile(mol, sdf_path)

pdbqt_path = "dataset/results/pdbqt/"
sdf_path = "dataset/results/sdf/"



for result in tqdm(os.listdir(pdbqt_path)):
    
    out_name = result.replace(".pdbqt", ".sdf")

    convert_pdbqt_to_sdf(pdbqt_path+result, sdf_path+out_name)