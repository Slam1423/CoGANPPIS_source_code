# CoGANPPIS_source_code
## Overview
This repository is the source code of CoGANPPIS, a coevolution-enhanced global attention neural network for protein-protein interaction site prediction. It contains the raw datasets, the feature generator, the pretrained model as well as the training scripts.

## Installation
```bash
git lfs clone git@github.com:Slam1423/CoGANPPIS_source_code.git
```

## Requirements:
- Linux
- python 3.6.9+
- biopython 1.79
- pybiolib 1.1.988
- pytorch 1.10.0
- scikit-learn 0.24.2
- numpy 1.19.5
- pandas 1.1.5
- scipy 1.5.4

## How to run
### Step1: Feature generation.

#### Hyperparameters:
`-f`: The dataset name.
usage example: `-f example`

`-d`: Protein database for Blast.
usage example: `-d nr`

`-n`: Number of sequences in multiple sequence alignments. (preferably >= 1000)
usage example: `-n 1000`

```bash
python3 generate_features.py -f dset422 -d nr -n 1000
```

### Step2: Training

```bash
python3 train.py
```

