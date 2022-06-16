# SMILES Mapper

## convert SMILES to map structure, or convert the map structure to SMILES

## Usage

```bash
 # in the parent folder of src
 python ./src/main.py
```
 + First, the **examples.csv** in folder "assets" will be read.
 + Then, **output.json** from decoder will be generated, containing the data of map structure like this:
  ```js
  // in output.json
    "7": [  // The 7th smiles in example.csv
        {
            "atom": "C",
            "index": 0, // THe first atom in smiles 
            "neighbors": [ // list of its neighbors
                {
                    "index": 1, 
                    "bond": 1, // bond number of the atom and its neighbors
                    "atom": "Mn"
                }
            ]
        },
        {
            "atom": "Mn",
            "index": 1,
            "neighbors": [
                {
                    "index": 0,
                    "bond": 1,
                    "atom": "C"
                },
                {
                    "index": 2,
                    "bond": 1,
                    "atom": "C"
                }
            ]
        },
        ...
  ```
 + Finally, **result.json** will be generated. This is the output of encoder, which take output.json as its input.
```js
  // in result.json
  {
    "1": "CCCCC(CC)OCC",
    "2": "CCOCc1ccccc1CC",
    "3": "CCOCc1cccc(c1)CC",
    "4": "CCOCCCNc1cc(C)ccc1N",
    "5": "c2ccc1c(c2)cccc1",
    "6": "C[13C]OCC[Fe]CC(SCC[Mn]C)CS[OH3+]",
    "7": "C[Mn]CC(C(C=S)[Fe]CCC=O)SC"
}
```

## Improvement

+ The signals **\ and @** are not already supported

## Announcement

+ **DO NOT** use this code in your paper, it is for learning and still under testing