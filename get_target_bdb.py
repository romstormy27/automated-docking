import xmltodict
import pprint
import json
import requests
import pubchempy as pcp
from tqdm import tqdm

"""

in order to get the target protein using bdb API,
we need to provide the API with SMILES from desired ligand/compound

"""

# get canonical smiles from pubchem using list of cid
def get_smiles(cid_list:list, save=False)->pd.DataFrame:

    cid_df = pd.DataFrame()

    for cid in cid_list:

        cid_df = pd.concat([cid_df, pcp.get_properties(["canonical_smiles"], cid, "cid", as_dataframe=True)])

    cid_df.reset_index(inplace=True)
    cid_df.columns = ["cid", "smiles"]

    if save:
        
        cid_df.to_csv("dataset/cid_smiles.csv", index=False)

    return cid_df


# use smiles to get target protein data from binding db
def get_target_data(cid_smiles:pd.DataFrame, similarity=0.9,save=False)->pd.DataFrame:

    df = pd.DataFrame()

    smiles_list = cid_smiles["smiles"].tolist()

    cid_list = cid_smiles["cid"].tolist()

    for smiles, cid in tqdm(zip(smiles_list, cid_list)):

        web_path = f"https://bindingdb.org/axis2/services/BDBService/getTargetByCompound?smiles={smiles}&cutoff={similarity}"

        bdb = requests.get(web_path)

        result = xmltodict.parse(bdb.content)

        try:

            tmp = result["bdb:getTargetByCompoundResponse"]["bdb:affinities"]

            for i in range(len(tmp)):

                tmp[i].update({"bdb:cid": str(cid)})

            tmp = pd.DataFrame(tmp)

            new_cols = [col.replace("bdb:", "") for col in tmp.columns]

            tmp.columns = new_cols

            df = pd.concat([df, tmp])

        except KeyError:

            print("targets for cid: ", cid, " are not found", end="")

    return df[["cid", "target", "species", "affinity", "affinity_type", "tanimoto"]]


if __name__=="__main__":

    from utils import load_config

    config = load_config()

    # load cid list from config files
    cid_list = config["ligand_cid_list"]

    # run get smiles function and store into dataframe
    cid_df = get_smiles(cid_list=cid_list, save=False)

    # run bdb api to get desired targets from a ligand and store into dataframe
    target_df = get_target_data(cid_smiles=cid_df)

        