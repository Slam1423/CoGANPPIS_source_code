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
- numpy 1.19.5
- pandas 1.1.5
- scipy 1.5.4

## How to run
First, you have to copy your dataset into `CoGANPPIS_package/raw_input_sequences/` with the name of `datasetname + _seq.txt` in the form of fasta. For example, your dataset name is `example`, then you should save the sequences into `example_seq.txt` in the form of fasta as follows:

```bash
>1acb_I
KSFPEVVGKTVDQAREYFTLHYPQYDVYFLPEGSPVTLDLRYNRVRVFYNPGTNVVNHVPHVG
>1ay7_A
DVSGTVCLSALPPEATDTLNLIASDGPFPYSQDGVVFQNRESVLPTQSYGYYHEYTVITPGARTRGTRRIITGEATQEDYYTGDHYATFSLIDQTC
>1c1y_B
SNTIRVFLPNKQRTVVNVRNGMSLHDCLMKALKVRGLQPECCAVFRLLHEHKGKKARLDWNTDAASLIGEELQVDFL

```

Now you can run the package with the following codes:

```bash
cd CoGANPPIS_package/
python3 main.py -f example -d nr -n 100
```

The predicted labels will be saved to `CoGANPPIS_package/predict_result.pkl` in the form of list in python, whose elements sequentially refer to the predicted labels of the residues and length is equal to the total number of residues of the input sequences. 

