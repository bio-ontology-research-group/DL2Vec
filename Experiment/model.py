import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import numpy as np




class Rank_model(nn.Module):
    def __init__(self, num_feature):
        super(Rank_model, self).__init__()

        self.model = nn.Sequential(
            nn.Linear(num_feature, 256),
            nn.Dropout(0.2),
            # nn.BatchNorm1d(256),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(256, 50),
            # nn.Dropout(0.5),
            # nn.BatchNorm1d(50),
            nn.LeakyReLU(0.2, inplace=True),
            # nn.Linear(256, 1),
            # nn.Sigmoid()
        )
        self.model2=nn.Sequential(
            nn.Linear(num_feature, 256),
            nn.Dropout(0.2),
            # nn.BatchNorm1d(256),
            nn.LeakyReLU(0.2, inplace=True),
            nn.Linear(256, 50),
            # nn.BatchNorm1d(50),
            # nn.Dropout(0.5),
            nn.LeakyReLU(0.2, inplace=True),
            # nn.Linear(256, 1),
            # nn.Sigmoid()
        )
        self.output_sig = nn.Sigmoid()

    def forward(self, gene_1,disease):

        g1 = self.model(gene_1).view(-1,1,50)
        d1 = self.model2(disease).view(-1,50,1)

        s1=torch.bmm(g1,d1)

        s1=s1.reshape(-1)



        out = self.output_sig(s1)

        return out

    def predict(self, data):


        gene=data[:,0,:]
        disease=data[:,1,:]


        gene_embed = self.model(gene)
        gene_embed=gene_embed.view(-1,1,50)


        disease_embed= self.model2(disease)
        disease_embed=disease_embed.view(-1,50,1)


        similarity=torch.bmm(gene_embed,disease_embed).view(-1)


        return similarity