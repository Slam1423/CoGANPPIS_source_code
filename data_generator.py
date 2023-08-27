import os
import time
import pickle
import torch as t
import numpy as np
from torch.utils import data
from config import DefaultConfig 


configs = DefaultConfig()


class dataSet(data.Dataset):
    def __init__(self,window_size,sequences_file=None,pssm_file=None, dssp_file=None, label_file=None, protein_list_file=None, train_MSA_file=None):
        super(dataSet,self).__init__()
        
        self.all_sequences = []
        for seq_file in sequences_file:
            with open(seq_file,"rb") as fp_seq:
               temp_seq  = pickle.load(fp_seq)
            self.all_sequences.extend(temp_seq)

        self.all_pssm = []
        for pm_file in pssm_file: 
            with open(pm_file,"rb") as fp_pssm:
                temp_pssm = pickle.load(fp_pssm)
            self.all_pssm.extend(temp_pssm)

        self.all_dssp = []
        for dp_file in dssp_file: 
            with open(dp_file,"rb") as fp_dssp:
                temp_dssp  = pickle.load(fp_dssp)
            self.all_dssp.extend(temp_dssp)

        self.all_label = []
        for lab_file in label_file: 
            with open(lab_file, "rb") as fp_label:
                temp_label = pickle.load(fp_label)
            self.all_label.extend(temp_label)

        with open(protein_list_file, "rb") as list_label:
            self.protein_list = pickle.load(list_label)
        
        
        self.all_dca = []
        self.MSA_features_matrix = []
        for dca_file in train_MSA_file: 
            with open(dca_file, "rb") as fp_dca:
                temp_dca = pickle.load(fp_dca)
            self.MSA_features_matrix.extend(temp_dca)

        self.max_seq_len = configs.max_sequence_length
        self.window_size = window_size


    def __getitem__(self,index):
        count,id_idx,ii,dset,protein_id,seq_length = self.protein_list[index]
        window_size = self.window_size
        id_idx = int(id_idx)
        win_start = ii - window_size
        win_end = ii + window_size
        seq_length = int(seq_length)
        label_idx = (win_start+win_end)//2
        msa_feature = np.array(self.MSA_features_matrix[id_idx][ii])
        
        idx = self.all_sequences[id_idx][ii]
        middle_data = []
        acid_one_hot = [0 for i in range(configs.seq_dim)]
        acid_one_hot[idx] = 1
        middle_data.extend(acid_one_hot)
        
        pssm_val = self.all_pssm[id_idx][ii]
        middle_data.extend(pssm_val)
        
        try:
            dssp_val = self.all_dssp[id_idx][ii]
            dssp_val.extend(0)
        except:
            dssp_val = [0 for i in range(configs.netsurf_dim)]
        middle_data.extend(dssp_val)
        middle_feature = np.stack(middle_data)
        
        
        all_seq_features = []
        seq_len = 0
        for idx in self.all_sequences[id_idx][:self.max_seq_len]:
            acid_one_hot = [0 for i in range(configs.seq_dim)]
            acid_one_hot[idx] = 1
            all_seq_features.append(acid_one_hot)
            seq_len += 1
        while seq_len<self.max_seq_len:
            acid_one_hot = [0 for i in range(configs.seq_dim)]
            all_seq_features.append(acid_one_hot)
            seq_len += 1

        all_pssm_features = self.all_pssm[id_idx][:self.max_seq_len]
        seq_len = len(all_pssm_features)
        while seq_len<self.max_seq_len:
            zero_vector = [0 for i in range(configs.pssm_dim)]
            all_pssm_features.append(zero_vector)
            seq_len += 1

        all_dssp_features = self.all_dssp[id_idx][:self.max_seq_len]
        seq_len = len(all_dssp_features)
        while seq_len<self.max_seq_len:
            zero_vector = [0 for i in range(configs.netsurf_dim)]
            all_dssp_features.append(zero_vector)
            seq_len += 1

        local_features = []
        labels = []
        while win_start<0:
            data = []
            acid_one_hot = [0 for i in range(configs.seq_dim)]
            data.extend(acid_one_hot)

            pssm_zero_vector = [0 for i in range(configs.pssm_dim)]
            data.extend(pssm_zero_vector)

            dssp_zero_vector = [0 for i in range(configs.netsurf_dim)]
            data.extend(dssp_zero_vector)

            local_features.extend(data)
            win_start += 1
       
        valid_end = min(win_end,seq_length-1)
        while win_start<=valid_end:
            data = []
            idx = self.all_sequences[id_idx][win_start]

            acid_one_hot = [0 for i in range(configs.seq_dim)]
            acid_one_hot[idx] = 1
            data.extend(acid_one_hot)


            pssm_val = self.all_pssm[id_idx][win_start]
            data.extend(pssm_val)

            try:
                dssp_val = self.all_dssp[id_idx][win_start]
            except:
                dssp_val = [0 for i in range(configs.netsurf_dim)]
            data.extend(dssp_val)

            local_features.extend(data)
            win_start += 1

        while win_start<=win_end:
            data = []
            acid_one_hot = [0 for i in range(configs.seq_dim)]
            data.extend(acid_one_hot)

            pssm_zero_vector = [0 for i in range(configs.pssm_dim)]
            data.extend(pssm_zero_vector)

            dssp_zero_vector = [0 for i in range(configs.netsurf_dim)]
            data.extend(dssp_zero_vector)

            local_features.extend(data)
            win_start += 1


        label = self.all_label[id_idx][label_idx]
        label = np.array(label,dtype=np.float32)

        all_seq_features = np.stack(all_seq_features)
        all_seq_features = all_seq_features[np.newaxis,:,:]
        all_pssm_features = np.stack(all_pssm_features)
        all_pssm_features = all_pssm_features[np.newaxis,:,:]

        all_dssp_features = np.stack(all_dssp_features)
        all_dssp_features = all_dssp_features[np.newaxis,:,:]
        local_features = np.stack(local_features)
        
        return all_seq_features,all_pssm_features,all_dssp_features,local_features,label, msa_feature, middle_feature
                

    def __len__(self):
    
        return len(self.protein_list)
