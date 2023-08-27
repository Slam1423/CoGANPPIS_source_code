import sys
from Bio.Blast import NCBIWWW
import biolib
import pickle
import numpy as np
import pandas as pd
import re
from plmDCA import plmDCA_main
import signal


class Timeout(Exception):
    pass

def handler(sig, frame):
    raise Timeout


dataset_name = sys.argv[2]
database = sys.argv[4]
hitlist_size = int(sys.argv[6])
input_file = 'dataset/' + dataset_name + '_seq.txt'

signal.signal(signal.SIGALRM, handler)

# raw sequence -> MSA files
file = open(input_file, 'r')
target_protein_file_list = file.readlines()
file.close()

blast_num = 20
blast_protein_num = 10

target_len_list = []
target_protein_list = []
target_protein_name_list = []
target_len_list1 = []
target_protein_list1 = []
target_protein_name_list1 = []
file = open(input_file, 'r')
target_file_list = file.readlines()
file.close()
for i in range(len(target_file_list)):
    if '>' in target_file_list[i]:
        target_protein_name_list1.append(target_file_list[i][1:-1])
        continue
    target_len_list1.append(len(target_file_list[i][:-1]))
    target_protein_list1.append(target_file_list[i][:-1])

fasta_list = []
for i in range(len(target_protein_file_list)):
    if blast_num * (i+1) < len(target_protein_file_list):
        fasta_list.append(''.join(target_protein_file_list[blast_num*i: blast_num*(i+1)]))
    else:
        fasta_list.append(''.join(target_protein_file_list[blast_num * i:]))
        break

total_cnt = 0
print('Start netsurfp...')
print('biolib.__version__: ' + biolib.__version__)
nsp3 = biolib.load('DTU/NetSurfP-3')
print('dataset/' + dataset_name + '_seq.txt')
nsp3_results = nsp3.cli(args='-i dataset/' + dataset_name + '_seq.txt')
nsp3_results.save_files("biolip_netsurfp/")
netsurf_df = pd.read_csv('biolip_netsurfp/results.csv', header=0)
asa_max = netsurf_df.iloc[:, 4].max()
sequence_data_list = []
faked_label_list = []
pssm_data_list = []
protein_netsurf_list = []
msa_feature_list = []
dset_list = []
print('Start Blastp...')
for iii in range(len(fasta_list)):
    fasta = fasta_list[iii]
    print('Blastp group: ' + str(iii) + ' (Mostly 10 sequences in a group)')
    flag = True
    while flag:
        signal.alarm(7200)
        try:
            result_handle = NCBIWWW.qblast('blastp', database, fasta, alignments=2000, hitlist_size=hitlist_size, format_type="Text")
            res_str = result_handle.read()
            if 'Query' in res_str and 'Sbjct' in res_str and 'ka-blk-sigma' in res_str:
                pass
            else:
                print('continue')
                continue
            file = open('msa_dir/protein_msa_multi_output' + str(iii) + '.txt', 'w')
            file.write(res_str)
            file.close()
            flag = False
        except Exception as e:
            print(e)
        signal.alarm(0)

    # MSA file -> pure MSA
    output_file = open('msa_dir/protein_msa_multi_output' + str(iii) + '.txt', 'r')
    raw_total_list = output_file.read()

    for i in range(len(target_protein_name_list1)):
        name = target_protein_name_list1[i]
        if 'Query= ' + name in raw_total_list:
            target_len_list.append(target_len_list1[i])
            target_protein_list.append(target_protein_list1[i])
            target_protein_name_list.append(name)

    total_list = raw_total_list.split('Query=')[1:]
    for ii in range(len(total_list)):
        each = total_list[ii]
        target_protein = target_protein_list[blast_protein_num * iii + ii]
        curList = each.split('Posted')[0].split('Score')
        length = target_len_list[blast_protein_num * iii + ii]
        msa_list = []
        whole_str_file = ''
        for i in range(len(curList)):
            cur_part = curList[i].split('\n')
            cur_msa_list = ['-' for u in range(length)]
            curcnt = 0
            for line in cur_part:
                if 'Query' in line:
                    lineList = re.split(r"[ ]+", line)
                    s = int(lineList[1])
                    t = int(lineList[-1])
                    curFatherSeq = lineList[2]
                if 'Sbjct' in line:
                    lineList = re.split(r"[ ]+", line)
                    cur_msa = lineList[2]
                    for j in range(len(cur_msa)):
                        if curFatherSeq[j] == '-':
                            continue
                        c_str = cur_msa[j]
                        cur_msa_list[curcnt] = c_str
                        curcnt += 1
            cur_msa_str = ''.join(cur_msa_list)
            msa_list.append(cur_msa_str)
            whole_str_file += cur_msa_str + '\n'
        file = open('msa_dir/protein_pure_msa_output' + str(ii + blast_protein_num * iii) + '.txt', 'w')
        file.write(whole_str_file)
        file.close()

# pure MSA -> DCA
print('Start computing DCA, PSSM...')
for ii in range(len(target_protein_name_list)):
    target_protein = target_protein_list[ii]
    plmDCA_main('msa_dir/protein_pure_msa_output' + str(ii) + '.txt', 'msa_dir/protein_pure_dca_output' + str(ii) + '.pkl', 0.1)
    output_file = open('msa_dir/protein_pure_dca_output' + str(ii) + '.pkl', 'rb')
    values = pickle.load(output_file)
    output_file.close()
    protein_seq = target_len_list[ii]
    dca_mat = np.ones((protein_seq, protein_seq))
    cnt = 0
    for i in range(protein_seq):
        for j in range(protein_seq):
            if i >= j:
                if i > j:
                    dca_mat[i][j] = 0
                continue
            dca_mat[i][j] = values[cnt]
            cnt += 1
    L = 500
    dca_list = []
    for i in range(protein_seq):
        if protein_seq > L:
            dca_list.append(dca_mat[i][:L].tolist())
        else:
            dca_list.append(dca_mat[i].tolist())
    protein_len = target_len_list[ii]
    if protein_len < L:
        for i in range(len(dca_list)):
            dca_list[i] += [0 for u in range(L - protein_len)]
    msa_feature_list.append(dca_list)

    # pure MSA -> PSSM
    print('pure MSA -> PSSM')
    pure_msa_file = open('msa_dir/protein_pure_msa_output' + str(ii) + '.txt', 'r')
    content = pure_msa_file.readlines()
    pure_msa_file.close()
    amino_acid_list = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    letter_to_idx = dict()
    for i in range(len(amino_acid_list)):
        letter_to_idx[amino_acid_list[i]] = i
    prime_seq = content[0][:-1]
    theLen = target_len_list[ii]
    matrix = [[0 for u in range(20)] for uu in range(theLen)]
    for i in range(len(content)):
        cur = content[i][:-1]
        for j in range(len(cur)):
            now = cur[j]
            if now not in letter_to_idx:
                continue
            idx = letter_to_idx[now]
            matrix[j][idx] += 1
    for i in range(theLen):
        for j in range(20):
            matrix[i][j] = int(np.log2(20 * (matrix[i][j] + 1) / (len(content) + 1)))
    pssm_data_list.append(matrix)

    # netsurfp
    cur_df = netsurf_df[netsurf_df['id'] == '>' + target_protein_name_list[ii]]
    cur_len = cur_df.shape[0]
    need_cols = [3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17]
    secondary_structure = cur_df.iloc[:, need_cols].values.tolist()
    for i in range(len(secondary_structure)):
        secondary_structure[i].append(target_len_list[ii] / 700)
        secondary_structure[i][1] /= asa_max
    protein_netsurf_list.append(secondary_structure)

    # sequence_data.pkl
    num_to_amino_acid = {
        0: 'ALA',
        1: 'CYS',
        2: 'ASP',
        3: 'GLU',
        4: 'PHE',
        5: 'GLY',
        6: 'HIS',
        7: 'ILE',
        8: 'LYS',
        9: 'LEU',
        10: 'MET',
        11: 'ASN',
        12: 'PRO',
        13: 'GLN',
        14: 'ARG',
        15: 'SER',
        16: 'THR',
        17: 'VAL',
        18: 'TRP',
        19: 'TYR'
    }

    amino_acid_to_letters = {
        'ALA': 'A',
        'ARG': "R",
        'ASN': 'N',
        'ASP': 'D',
        'CYS': 'C',
        'GLN': 'Q',
        'GLU': 'E',
        'GLY': 'G',
        'HIS': 'H',
        'ILE': 'I',
        'LEU': 'L',
        'LYS': 'K',
        'MET': 'M',
        'PHE': 'F',
        'PRO': 'P',
        'SER': 'S',
        'THR': 'T',
        'TRP': 'W',
        'TYR': 'Y',
        'VAL': 'V'
    }

    letter_to_number = dict()
    for key in num_to_amino_acid:
        mapped_letter = amino_acid_to_letters[num_to_amino_acid[key]]
        letter_to_number[mapped_letter] = key
    seq_list = []
    fake_labels = []
    for i in range(theLen):
        fake_labels.append(0)
        if target_protein[i] in letter_to_number:
            seq_list.append(letter_to_number[target_protein[i]])
        else:
            seq_list.append(letter_to_number['N'])
    sequence_data_list.append(seq_list)
    faked_label_list.append(fake_labels)

    # dset_list.pkl
    for i in range(theLen):
        cur = (total_cnt, ii, i, dataset_name, target_protein_name_list[ii], str(theLen))
        dset_list.append(cur)
        total_cnt += 1

pkl_file = open('data_cache/' + dataset_name + '_sequence_data.pkl', 'wb')
pickle.dump(sequence_data_list, pkl_file)
pkl_file.close()
pkl_file = open('data_cache/' + dataset_name + '_pssm_data.pkl', 'wb')
pickle.dump(pssm_data_list, pkl_file)
pkl_file.close()
pkl_file = open('data_cache/' + dataset_name + '_netsurf_ss_14.pkl', 'wb')
pickle.dump(protein_netsurf_list, pkl_file)
pkl_file.close()
pkl_file = open('data_cache/' + dataset_name + '_MSA_features_1.pkl', 'wb')
pickle.dump(msa_feature_list, pkl_file)
pkl_file.close()
pkl_file = open('data_cache/' + dataset_name + '_list.pkl', 'wb')
pickle.dump(dset_list, pkl_file)
pkl_file.close()
