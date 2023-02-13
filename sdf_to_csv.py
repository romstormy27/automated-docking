import pandas as pd
from rdkit import Chem
import argparse
from collections import defaultdict


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-input", type=str, required=True)
    parser.add_argument("-output", type=str, required=True)
    parser.add_argument("-save_name", action='store_true', help="store header line as _Name")
    args = parser.parse_args()

    # SDFの読み込み(1回目はパラメータ名を全て読み取る)
    sdf_sup = Chem.SDMolSupplier(args.input)
    Props = []
    if args.save_name:
        Props.append("_Name")

    for mol in sdf_sup:
        for name in mol.GetPropNames():
            if name not in Props:
                Props.append(name)

    # データを格納するディクショナリ
    param_dict = defaultdict(list)

    # SDFの読み込み(2回目は化合物のパラメータを取得。なければエラー)
    sdf_sup = Chem.SDMolSupplier(args.input)
    for mol in sdf_sup:
        # 名前の取得
        for name in Props:
            if mol.HasProp(name):
                param_dict[name].append(mol.GetProp(name))
            else:
                param_dict[name].append(None)

    #　pandasで一気に変換
    df = pd.DataFrame(data=param_dict)
    df.to_csv(args.output, index=False)

if __name__ == "__main__":
    main()