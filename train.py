import os
import time
import pickle
import numpy as np
import torch
import random
from torch import nn
import torch.utils.data.sampler as sampler
from config import DefaultConfig
from ppi_model import PPIModel
import data_generator
from evaluation import compute_roc, compute_aupr, compute_mcc, micro_score, acc_score, compute_performance
from predict import make_prediction


configs = DefaultConfig()
global train_data
global model_save_path


def train_epoch(model, loader, optimizer, epoch, all_epochs, print_freq=100):

    model.train()

    end = time.time()
    for batch_idx, (seq_data, pssm_data, dssp_data, local_data, label, msa_file, middle_fea) in enumerate(loader):
        end = time.time()
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
        output = output.view(shapes[0] * shapes[1])
        for_weights = target_var.float()
        weight_var = for_weights + configs.weighted_loss
        loss = torch.nn.functional.binary_cross_entropy(output, target_var, weight_var).cuda()

        batch_size = label.size(0)

        # compute gradient and do SGD step
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        # print stats
        if batch_idx % print_freq == 0:
            res = '\t'.join([
                'Epoch: [%d/%d]' % (epoch + 1, all_epochs),
                'Iter: [%d/%d]' % (batch_idx + 1, len(loader)),
                'Time %.3f' % (time.time() - end),
                'Loss %.4f' % (loss.item())])
            print(res)

    return time.time() - end, loss.item()


def eval_epoch(model, loader, print_freq=10, is_test=True, test_file=None):

    # Model on eval mode
    model.eval()

    all_trues = []
    all_preds = []
    all_gos = []
    end = time.time()
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

        # compute output
        output = model(seq_var, dssp_var, pssm_var, local_var, msa_var, middle_var)
        shapes = output.data.shape
        output = output.view(shapes[0] * shapes[1])

        for_weights = target_var.float()

        weight_var = for_weights + configs.weighted_loss
        loss = torch.nn.functional.binary_cross_entropy(output, target_var, weight_var).cuda()

        # measure accuracy and record loss
        batch_size = label.size(0)
        # print stats
        batch_time = time.time() - end
        all_trues.append(label.numpy())
        all_preds.append(output.data.cpu().numpy())

    all_trues = np.concatenate(all_trues, axis=0)
    all_preds = np.concatenate(all_preds, axis=0)

    auc = compute_roc(all_preds, all_trues)
    aupr = compute_aupr(all_preds, all_trues)
    f_max, p_max, r_max, t_max, predictions_max, threshold_max = compute_performance(all_preds, all_trues)
    acc_val = acc_score(predictions_max, all_trues)
    mcc = compute_mcc(predictions_max, all_trues)
    
    return batch_time, loss.item(), acc_val, f_max, p_max, r_max, auc, aupr, t_max, mcc, threshold_max


def train(dataset_list, model, train_data_set, save, n_epochs=30, batch_size=256, learning_rate=0.0001, train_file=None, valid_file=None, test_file=None):

    dataset_str = '_'.join(dataset_list)

    with open(train_file, "rb") as fp:
        train_list = pickle.load(fp)

    with open(valid_file, "rb") as fp_val:
        valid_list = pickle.load(fp_val)

    with open(test_file, "rb") as fp_val:
        test_list = pickle.load(fp_val)

    train_index = train_list
    eval_index = valid_list
    test_index = test_list
    train_samples = sampler.SubsetRandomSampler(train_index)
    eval_samples = sampler.SubsetRandomSampler(eval_index)
    test_samples = sampler.SubsetRandomSampler(test_index)

    train_loader = torch.utils.data.DataLoader(train_data_set, batch_size=batch_size, sampler=train_samples, pin_memory=True, num_workers=4, drop_last=False)
    valid_loader = torch.utils.data.DataLoader(train_data_set, batch_size=batch_size, sampler=eval_samples, pin_memory=True, num_workers=4, drop_last=False)
    test_loader = torch.utils.data.DataLoader(train_data_set, batch_size=batch_size, sampler=test_samples, pin_memory=True, num_workers=4, drop_last=False)

    if torch.cuda.is_available():
        model = model.cuda()

    model_wrapper = model

    optimizer = torch.optim.Adam(filter(lambda p: p.requires_grad, model_wrapper.parameters()), lr=learning_rate)

    # Train model
    best_F = 0
    threadhold = 0
    count = 0
    for epoch in range(n_epochs):
        _, train_loss = train_epoch(
            model=model_wrapper,
            loader=train_loader,
            optimizer=optimizer,
            epoch=epoch,
            all_epochs=n_epochs,
        )
        _, valid_loss, acc, f_max, p_max, r_max, auc, aupr, t_max, mcc_max, threshold_max = eval_epoch(
            model=model_wrapper,
            loader=valid_loader,
            is_test=(not valid_loader),
            test_file=test_file
        )

        print('Validation:\nvalid_loss:%0.5f,threadhold:%0.6f' % (valid_loss, t_max))
        model_save_path = os.path.join(save, dataset_str + '_training_epoch_' + str(epoch) + '.dat')
        torch.save(model.state_dict(), model_save_path)
        
        predicted_f = make_prediction(model_save_path, train_data, threshold_max)
        if predicted_f > best_F:
            count = 0
            best_F = predicted_f
            THREADHOLD = t_max
            print("new best F_value:{0}(threadhold:{1})".format(predicted_f, THREADHOLD))
        else:
            count += 1
            if count >= 5:
                return None


def demo(train_data, save=None, window_size=3, epochs=30, pretrained_result=None):
    train_sequences_file = ['data_cache/{0}_sequence_data.pkl'.format(key) for key in train_data]
    train_dssp_file = ['data_cache/{0}_netsurf_ss_14.pkl'.format(key) for key in train_data]
    train_pssm_file = ['data_cache/{0}_pssm_data.pkl'.format(key) for key in train_data]
    train_label_file = ['data_cache/{0}_label.pkl'.format(key) for key in train_data]
    train_MSA_file = ['data_cache/{0}_MSA_features_1.pkl'.format(key) for key in train_data]
    
    if train_data == ['dset422']:
        all_list_file = 'data_cache/dset422_list.pkl'
        train_list_file = 'data_cache/dset422_train_list.pkl'
        valid_list_file = 'data_cache/dset422_val_list.pkl'
        test_list_file = 'data_cache/dset422_test_list.pkl'
    else:
        all_list_file = 'data_cache/dset448_list.pkl'
        train_list_file = 'data_cache/dset448_train_list.pkl'
        valid_list_file = 'data_cache/dset448_val_list.pkl'
        test_list_file = 'data_cache/dset448_test_list.pkl'

    batch_size = configs.batch_size

    print('batch_size:' + str(batch_size))

    train_dataSet = data_generator.dataSet(window_size, train_sequences_file, train_pssm_file, train_dssp_file, train_label_file, all_list_file, train_MSA_file)
    
    class_nums = 1
    model = PPIModel(class_nums, window_size)
    model_file = "pretrained_model.dat"
    model.load_state_dict(torch.load(model_file))

    for param in model.named_parameters():
        if ('DNN' not in param[0]) and ('outLayer' not in param[0]):
            param[1].requires_grad = False
        else:
            param[1].requires_grad = True

    train(train_data, model=model, train_data_set=train_dataSet, save=save,
          n_epochs=epochs, batch_size=batch_size, 
          train_file=train_list_file, valid_file=valid_list_file, test_file=test_list_file)
    print('Done!')


if __name__ == '__main__':

    path_dir = "./saved_models/ppi_model"
    train_data = ['dset422']
    
    if not os.path.exists(path_dir):
        os.makedirs(path_dir)

    demo(train_data, path_dir)
    