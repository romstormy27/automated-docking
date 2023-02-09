from vina import Vina
import pandas as pd
from utils import load_config
import sys
import subprocess

def get_interaction_map(path:str)->(dict):

    # read data
    df = pd.read_csv(path)

    # convert into dict and group by receptor
    tmp = df.reset_index()
    dock_dict = {rec: lig.tolist() for rec, lig in tmp.groupby("uniprot")["cid"]}
    dock_dict

    return dock_dict

def run_vina_docking(interaction_map:dict, receptor_path:str, ligand_path:str, output_path:str):

    for receptor in interaction_map:
        for ligand in interaction_map[receptor]:

            # print(receptor)
            # print(ligand)

            v = Vina(sf_name='vina', seed=42, verbosity=1)

            print("Performing docking simulation...")
            # subprocess.run(["echo", "-e", "Performing docking simulation...\n"])

            print()
            print(f"RECEPTOR: {receptor} <---> LIGAND: {ligand}")
            # subprocess.run(["echo", "-e", f"RECEPTOR: {receptor} <---> LIGAND: {ligand}\n"])
            print()

            v.set_receptor(rigid_pdbqt_filename=receptor_path+str(receptor)+".pdbqt")
            v.set_ligand_from_file(ligand_path+str(ligand)+".pdbqt")

            v.compute_vina_maps(center=[0,0,0], box_size=[80, 80, 80])

            # Score the current pose
            energy = v.score()
            print(f'Score before minimization: {energy[0]} (kcal/mol)')

            # Minimized locally the current pose
            energy_minimized = v.optimize()
            print(f'Score after minimization : {energy_minimized[0]} (kcal/mol)')
            # v.write_pose('1ao6_ligand_minimized.pdbqt', overwrite=True)

            # Dock the ligand
            v.dock(exhaustiveness=32, n_poses=20)
            v.write_poses(output_path+str(receptor)+"_"+str(ligand)+".pdbqt", n_poses=5, overwrite=True)

            # subprocess.run(["echo", "-e", "\n"])
            # subprocess.run(["echo", "=============================================================="])
            print("==============================================================")

if __name__=="__main__":

    config = load_config()

    # load all necessary path

    interaction_data_path = config["interaction_data_path"]

    receptor_data_path = config["receptor_pdbqt_path"]

    ligand_data_path = config["ligand_pdbqt_path"]

    docking_output_path = config["docking_output_path"]

    # get interaction map
    dock_dict = get_interaction_map(interaction_data_path)

    # run docking
    run_vina_docking(interaction_map=dock_dict, receptor_path=receptor_data_path,
                     ligand_path=ligand_data_path, output_path=docking_output_path)