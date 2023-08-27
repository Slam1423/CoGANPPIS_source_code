import os
import time
import sys
import torch as t
from torch import nn
import torch.nn.functional as F
from torch.autograd import Variable
import pickle
from config import DefaultConfig


sys.path.append("../")
configs = DefaultConfig()


class global_attention(nn.Module):
    def __init__(self):
        super().__init__()
        self.q1 = nn.Parameter(t.rand(2*(configs.pssm_dim + configs.seq_dim + configs.netsurf_dim) + 1))
        self.linear1 = nn.Linear((configs.pssm_dim + configs.seq_dim + configs.netsurf_dim), (configs.pssm_dim + configs.seq_dim + configs.netsurf_dim))
        self.linear2 = nn.Linear((configs.pssm_dim + configs.seq_dim + configs.netsurf_dim), (configs.pssm_dim + configs.seq_dim + configs.netsurf_dim))
        self.layer_norm = nn.LayerNorm(configs.max_sequence_length)
        self.LeakyRelu=nn.LeakyReLU(negative_slope=0.01)
        self.softmax=nn.Softmax(dim=-1)

    def forward(self, features, mid_feature, msa_features):
        features = features.squeeze(1)
        s = t.mean(features, dim=1).unsqueeze(1).repeat(1, configs.max_sequence_length, 1)
        target_s = mid_feature.unsqueeze(1).repeat(1, configs.max_sequence_length, 1)
        
        y = t.cat((self.linear1(t.mul(s, features)), self.linear2(t.mul(target_s, features))), dim=-1) 
        y = t.cat((y, msa_features.unsqueeze(-1)), dim=-1)
        y = self.LeakyRelu(y)
        x = y * self.q1
        x = x.sum(dim=-1)
        x = self.layer_norm(x)
        attention = self.softmax(x).unsqueeze(-1).repeat(1, 1, configs.pssm_dim + configs.seq_dim + configs.netsurf_dim)
        output = t.mul(attention, features).sum(dim=1)
        
        return output


class PPIModel(nn.Module):
    def __init__(self,class_nums,window_size):
        super(PPIModel,self).__init__()
        global configs

        seq_dim = configs.seq_dim*configs.max_sequence_length
        
        self.seq_layers = nn.Sequential()
        self.seq_layers.add_module("seq_embedding_layer",
        nn.Linear(seq_dim,seq_dim))
        self.seq_layers.add_module("seq_embedding_ReLU",
        nn.ReLU())
        
        self.atten = global_attention()

        seq_dim = configs.seq_dim
        netsurf_dim = configs.netsurf_dim
        pssm_dim = configs.pssm_dim
        input_dim = 464
        self.bn0 = nn.BatchNorm1d(1)

        self.DNN1 = nn.Sequential()
        self.DNN1.add_module("DNN_layer1",
                            nn.Linear(input_dim, 256))
        self.DNN1.add_module("ReLU1",
                            nn.ReLU())
        
        self.DNN2 = nn.Sequential()
        self.DNN2.add_module("DNN_layer2",
                            nn.Linear(256, 64))
        self.DNN2.add_module("ReLU2",
                            nn.ReLU())
                            
        self.DNN3 = nn.Sequential()
        self.DNN3.add_module("DNN_layer3",
                            nn.Linear(64,8))
        self.DNN3.add_module("ReLU2",
                            nn.ReLU())                    
        

        self.outLayer = nn.Sequential(
            nn.Linear(8, class_nums),
            nn.Sigmoid())
            
            
        in_channel = 1
        hidden_channels = 1
        self.conv1 = nn.Sequential()
        self.conv1.add_module("conv1",
                              nn.Conv1d(in_channel, hidden_channels,
                                        kernel_size=13, stride=13)
                              )
        self.conv1.add_module("ReLU", nn.ReLU())
        self.conv1.add_module("pooling1", nn.MaxPool1d(kernel_size=3, stride=3))
        
        self.conv2 = nn.Sequential()
        self.conv2.add_module("conv2",
                              nn.Conv1d(in_channel, hidden_channels,
                                        kernel_size=15, stride=15)
                              )
        self.conv2.add_module("ReLU", nn.ReLU())
        self.conv2.add_module("pooling2", nn.MaxPool1d(kernel_size=3, stride=3))
        
        self.conv3 = nn.Sequential()
        self.conv3.add_module("conv3",
                              nn.Conv1d(in_channel, hidden_channels,
                                        kernel_size=17, stride=17)
                              )
        self.conv3.add_module("ReLU", nn.ReLU())
        self.conv3.add_module("pooling1", nn.MaxPool1d(kernel_size=3, stride=3))
        
        self.dropout_1 = nn.Dropout(p=configs.dropout_rate)
        self.dropout_2 = nn.Dropout(p=configs.dropout_rate)
        self.dropout_3 = nn.Dropout(p=configs.dropout_rate)
        
        self.W_a = nn.Parameter(t.rand(98, 98))
        self.v_a = nn.Parameter(t.rand(98))
        
    
    def forward(self,seq,netsurf,pssm,local_features, msa_features, middle_features):
        
        
        shapes = seq.data.shape
        features = seq.view(shapes[0],shapes[1]*shapes[2]*shapes[3])
        features = self.seq_layers(features)
        features = features.view(shapes[0],shapes[1],shapes[2],shapes[3])
        features = t.cat((features,netsurf,pssm),3)
        features = self.atten(features, middle_features, msa_features)
        msa_features = msa_features.unsqueeze(1)
        msa_features = self.bn0(msa_features)
        msa_features1 = self.conv1(msa_features).squeeze(1)
        msa_features2 = self.conv2(msa_features).squeeze(1)
        msa_features3 = self.conv3(msa_features).squeeze(1)
        msa_features = t.cat((msa_features1, msa_features2, msa_features3), dim=-1)
        features = t.cat((features, local_features, msa_features), 1)
        features = self.dropout_1(self.DNN1(features))
        features = self.dropout_2(self.DNN2(features))
        features = self.dropout_3(self.DNN3(features))
        features = self.outLayer(features)

        return features
    
    
    def load(self,path):
    
        self.load_state_dict(t.load(path))
        
    def save(self,name=None):

        
        if name is None:
            prefix = ""
            name = time.strftime("%y%m%d_%H:%M:%S.pth".format(prefix))
            
        t.save(self.state_dict(),name)
        return name

