import os
import time
import pickle
import numpy as np
import torch
from torch.optim import lr_scheduler
from torch.nn.init import xavier_normal,xavier_normal_
from torch import nn
import torch.utils.data.sampler as sampler
from config import DefaultConfig
from ppi_model import PPIModel
import data_generator
from evaluation import compute_roc, compute_aupr, compute_mcc, micro_score,acc_score, compute_performance


configs = DefaultConfig()


def test(model, loader,path_dir,threshold):

    # Model on eval mode
    model.eval()
    length = len(loader)
    result = []
    all_trues = []


    for batch_idx, (seq_data, pssm_data, dssp_data, local_data, label, msa_file, middle_fea) in enumerate(loader):
    
        # Create vaiables
        with torch.no_grad():
            if torch.cuda.is_available():
                seq_var = torch.autograd.Variable(seq_data.cuda().float())
                pssm_var = torch.autograd.Variable(pssm_data.cuda().float())
                dssp_var = torch.autograd.Variable(dssp_data.cuda().float())
                local_var = torch.autograd.Variable(local_data.cuda().float())
                target_var = torch.autograd.Variable(label.cuda().float())
                msa_var = torch.autograd.Variable(msa_file.cuda().float())
                middle_var = torch.autograd.Variable(middle_fea.cuda().float())
            else:
                seq_var = torch.autograd.Variable(seq_data.float())
                pssm_var = torch.autograd.Variable(pssm_data.float())
                dssp_var = torch.autograd.Variable(dssp_data.float())
                local_var = torch.autograd.Variable(local_data.float())
                target_var = torch.autograd.Variable(label.float())
                msa_var = torch.autograd.Variable(msa_file.float())
                middle_var = torch.autograd.Variable(middle_fea.float())

        output = model(seq_var, dssp_var, pssm_var, local_var, msa_var, middle_var)
        shapes = output.data.shape
        output = output.view(shapes[0]*shapes[1])
        result.append(output.data.cpu().numpy())
        all_trues.append(label.numpy())

    #caculate
    all_trues = np.concatenate(all_trues, axis=0)
    all_preds = np.concatenate(result, axis=0)

    auc = compute_roc(all_preds, all_trues)
    aupr = compute_aupr(all_preds, all_trues)
    f_max, p_max, r_max, t_max, predictions_max, threshold_for_predict = compute_performance(all_preds,all_trues, True, threshold)
    acc = acc_score(predictions_max,all_trues)
    mcc = compute_mcc(predictions_max, all_trues)

    print( 'Test:\nacc:%0.6f,F_value:%0.6f, precision:%0.6f,recall:%0.6f,auc:%0.6f,aupr:%0.6f,mcc:%0.6f,threadhold:%0.6f\n' % (
        acc, f_max, p_max, r_max,auc, aupr,mcc,threshold_for_predict))

    predict_result = {}
    predict_result["pred"] = all_preds
    predict_result["label"] = all_trues
    result_file = "{0}/predict_result.pkl".format(path_dir)

    with open(result_file,"wb") as fp:
        pickle.dump(predict_result,fp)
    return f_max


def predict(model_file,test_data,window_size,path_dir,threshold):
    test_sequences_file = ['data_cache/{0}_sequence_data.pkl'.format(key) for key in test_data]
    test_dssp_file = ['data_cache/{0}_netsurf_ss_14.pkl'.format(key) for key in test_data]
    test_pssm_file = ['data_cache/{0}_pssm_data.pkl'.format(key) for key in test_data]
    test_label_file = ['data_cache/{0}_label.pkl'.format(key) for key in test_data]
    test_MSA_file = ['data_cache/{0}_MSA_features_1.pkl'.format(key) for key in test_data]
    
    if test_data == ['dset422']:
        test_list_file = 'data_cache/dset422_test_list.pkl'
        all_list_file = 'data_cache/dset422_list.pkl'
    else:
        test_list_file = 'data_cache/dset448_test_list.pkl'
        all_list_file = 'data_cache/dset448_list.pkl'
    
    # parameters
    batch_size = configs.batch_size

    # Datasets
    test_dataSet = data_generator.dataSet(window_size, test_sequences_file, test_pssm_file, test_dssp_file, test_label_file,
                                             all_list_file, test_MSA_file)
    # Models
    
    with open(test_list_file,"rb") as fp:
        test_list = pickle.load(fp)

    test_samples = sampler.SubsetRandomSampler(test_list)
    test_loader = torch.utils.data.DataLoader(test_dataSet, batch_size=batch_size,
                                              sampler=test_samples, pin_memory=False,
                                               num_workers=0,drop_last=False)

    # Models
    class_nums = 1
    model = PPIModel(class_nums,window_size)
    model.load_state_dict(torch.load(model_file))
    model = model.cuda()
    return test(model, test_loader,path_dir,threshold)


def make_prediction(model_file_name, dataset_list, threshold):

    window_size = 3
    path_dir = "./saved_models/ppi_model"
    datas = dataset_list

    if not os.path.exists(path_dir):
        os.makedirs(path_dir)
    
    model_file = model_file_name
    return predict(model_file,datas,window_size,path_dir,threshold)

