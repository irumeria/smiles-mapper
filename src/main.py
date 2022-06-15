import numpy as np
import json
import pandas as pd
import re
from decoder import *
from encoder import *


def load_smiles_from_csv(data_file, max_len):
    df = pd.read_csv(data_file)
    df.iloc[:, 0] = df.iloc[:, 0].str.strip()
    df = df[df.iloc[:, 0].str.len() <= max_len]
    smiles = df.iloc[:, 0].tolist()
    return smiles

def test_encoder():
    with open('output.json', 'r+', encoding='UTF-8') as f:
        nodeMaps = json.load(f)
    nodeMaps = nodeMaps.items()
    smilesJSON = {}
    counter = 0
    for key, nodeMap in nodeMaps:
        counter += 1
        print("turn",key)
        smilesJSON[counter] = map_to_smiles(nodeMap)
        
    with open('result.json', 'w',encoding='UTF-8') as f:
        f.write(json.dumps(smilesJSON, ensure_ascii=False))

def test_decoder():
    smiles = np.array(load_smiles_from_csv("./assets/examples.csv", 120))
    print(smiles)
    nodeMaps = {}
    counter = 0
    for sm in smiles:
        counter += 1
        print(sm)
        nodeMaps[str(counter)] = smiles_to_map(sm)
    with open('output.json', 'w',encoding='UTF-8') as f:
        f.write(json.dumps(nodeMaps, ensure_ascii=False))

if __name__ == "__main__":
    test_decoder()
    test_encoder()