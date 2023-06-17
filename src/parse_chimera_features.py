import os
import re
import glob 
import argparse
from typing import Dict
from pathlib import Path

import numpy as np
import pandas as pd

def create_arg_parse():
    # Create the parser
    parser = argparse.ArgumentParser(
        description='Concat the features extracted by chimera.')

    # Add the arguments
    parser.add_argument('-ip', '--input_path',
                        type=str,
                        required=True,
                        help='Path to the folder with contains the txt files to concat')

    # Add the arguments
    parser.add_argument('-op', '--output_path',
                        type=str,
                        help='Path to save concataned csv file.')

    args = parser.parse_args()

    return args.input_path, args.output_path

def create_data(protein_path:str):

    data = {}

    input_path = os.path.join(protein_path, "*.txt")

    for file_path in glob.glob(input_path):
        
        with open(file_path, "r", encoding="utf-8") as f:

            # skipping first two lines, which are the header
            conteudo = f.readlines()[2:]

            # processing the input to get only the relevant information
            # which is [[number, atribute]]
            conteudo = np.array(["-".join([re.sub('[:\nA-B]', '', 
                                                x) for x in y.split("\t")[1:]]).split(".-") for y in conteudo])

            
            file_name, _= os.path.splitext(os.path.basename(file_path))
            
            data[file_name] = conteudo

    return data

def create_table(data:Dict[str,np.array]):

    # chosing first atribute to start table to be merged with others
    df = pd.DataFrame({"number" : data["kdhydro"][:, 0],
                    "kdhydro": data["kdhydro"][:, 1]}, 
                    dtype=np.number)

    for k, v in data.items():
        if k not in df.columns:
            df = df.merge(pd.DataFrame({"number" : v[:, 0],
                                        k: v[:, 1]}, 
                                    dtype=np.number),
                        on="number", how="outer")
            
    df["number"] = df["number"].astype(int)

    return df

def create_save_location(path:str, protein_name:str):

    if not path:
        Path("results").mkdir(parents=True, exist_ok=True)
        output_path = os.path.join("results", f"{protein_name}_structure.csv")
    
    elif os.path.isdir(path):
        output_path = os.path.join(path, f"{protein_name}_structure.csv")

    else:
        output_path = path

    return output_path

def main():

    protein_path, output_path = create_arg_parse()

    protein_name = os.path.basename(os.path.realpath(protein_path))

    data = create_data(protein_path)
    df = create_table(data)

    save_path = create_save_location(output_path, protein_name)
    print(f"Saved in: {save_path}")

    df.to_csv(save_path, index=False)



if __name__ == "__main__":
    main()