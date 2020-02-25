from __future__ import print_function
import json
import numpy as np
import pickle as pkl
import gensim
from networkx.readwrite import json_graph
from argparse import ArgumentParser

from model import Rank_model
import torch
from torch._C import *
import random
from torch.utils.data import DataLoader
from torch.utils.data import TensorDataset
from torch.autograd import Variable
import torch.optim as optim
from scipy.stats import rankdata
import sys
import random
random.seed(42)
import argparse
import sys
import os


parser = argparse.ArgumentParser(description = "pointwise ranking model to prioritize the candidates")

parser.add_argument("-associations", nargs = '?', metavar = "association files ", type = str,
                     help = "File contains the associations which are in dictionary format",default='')
parser.add_argument("-ranking_candidate", nargs ="?", metavar ="ranking candidates", type=str,
                        help = "files contains the candidates that need to be prioritized, and it is in a list format")
parser.add_argument("-embedding_model", nargs="?", metavar ="the embedding models", default="",
                    help = " files contain the model that can generate the embedding for each entity")

args = parser.parse_args()

associations = args.associations
ranking_candidate = args.ranking_candidate
embedding_model = args.embedding_model

if (associations is '' ):
	print ("\nError:Mandatory associations file that contains the association with the candidates, For help, run: python pointwise_ranking.py --help\n")
	sys.exit()
if (ranking_candidate is ''):
	print ("\nError: candidate file is missing. For help, run: python pointwise_ranking.py --help\n")
	sys.exit()
if (embedding_model is ''):
	print ("\nError: embedding model file is missing. For help, run: python pointwise_ranking.py --help\n")
	sys.exit()





negative_number=20


def negative_sampling(disease_gene,gene_set,central_disease):


    positive_genes=disease_gene[central_disease]
    negative_genes=[]
    while(len(negative_genes)<negative_number):
        gene=random.choices(gene_set)
        if gene not in positive_genes:
            negative_genes.append(gene[0])

    return negative_genes


def generate_label(positive_g,negative_genes):
    label=[]

    train_gene=[]
    negative_gene=[]
    train_gene.append(positive_g)
    label.append(1)
    for gene in negative_genes:


        train_gene.append(gene)
        label.append(0)


        # if label_select==1:
        #
        #     positive_gene.append(positive_g)
        #     negative_gene.append(gene)
        # else:
        #
        #     negative_gene.append(gene)
        #     positive_gene.append(positive_g)

    return train_gene,label


def generate_embedding(embed_dic,entities):

    embed_list=[]
    for entity in entities:

        # id=id_mapping[entity]

        embed_list.append(embed_dic[entity])
    return embed_list


def convert_to_torch_format(data):
    data=np.array(data,dtype="float64")
    return torch.from_numpy(data).double()


def generate_train_data(embed_dic,disease_gene,gene_set,disease_set):
    g1 = []

    d = []
    y = []

    for disease in disease_set:
        genes = disease_gene[disease]
        for gene in genes:
            negative_genes = negative_sampling(disease_gene, gene_set, disease)
            train_gene, label = generate_label(gene, negative_genes)

            train_gene = generate_embedding(embed_dic, train_gene)


            disease_embedding = generate_embedding(embed_dic, [disease]*(negative_number+1))

            g1.extend(train_gene)

            d.extend(disease_embedding)
            y.extend(label)
    g1 = convert_to_torch_format(g1)

    d = convert_to_torch_format(d)
    y = convert_to_torch_format(y)
    return g1, d,y


def load_data(embed_dic,disease_gene,gene_set,test_diseases):

    disease_set=[key for key in disease_gene.keys()]
    #train_diseases=disease_set[:int(len(disease_genes)*0.8)]
    train_diseases=[]
    for disease in disease_gene.keys():
        if disease in test_diseases:

            pass

        else:
            train_diseases.append(disease)
    g1,d,y=generate_train_data(embed_dic,disease_gene,gene_set,train_diseases)




    return g1,d,y,test_diseases,train_diseases


def train(model_,opt_,criterion_):
    model_.train()
    loader = TensorDataset(g1, d, y)
    loader_dataset = DataLoader(loader, batch_size=125, shuffle=True)
    losses = 0

    for gene1_,disease_, y_label in loader_dataset:
        gene1_batch = Variable(gene1_)

        disease_batch = Variable(disease_)

        predict_y = model_(gene1_batch,disease_batch)
        loss = criterion_(predict_y, y_label)

        opt_.zero_grad()
        loss.backward()
        opt_.step()
        losses += loss.item()

    return model, losses

def look_up_embed(embed_dic_,entity):

    entity_emb=embed_dic_[entity]

    return entity_emb

def top_k(evaluation,k):
    scores=[]
    for key in evaluation.keys():
        scores.extend(evaluation[key])

    top_k_result=[]
    for data in scores:
        if data <k:
            top_k_result.append(data)
    rank=len(top_k_result)/len(scores)
    return rank



def evaluation(model_,evaluation_data,embed_dic_,gene_set_):
    model_.eval()
    validation_rank = dict()
    for disease in evaluation_data:
        positive_genes=disease_genes[disease]

        disease_vec=look_up_embed(embed_dic_,disease)


        testing_list = []
        for gene in gene_set_:
            gene_vec = look_up_embed(embed_dic_,gene)
            testing_list.append([gene_vec,disease_vec])

        validation_list = []
        for gene in positive_genes:
            gene_vec = look_up_embed(embed_dic_,gene)
            validation_list.append([gene_vec,disease_vec])

        testing_list=convert_to_torch_format(testing_list)
        validation_list=convert_to_torch_format(validation_list)
        testing_list=Variable(testing_list)
        validation_list=Variable(validation_list)

        test_result=sorted(model_.predict(testing_list).data,reverse=True)
        valid_result=model_.predict(validation_list).data

        rank_list = []


        for validation in valid_result:

            for index in range(len(test_result)):

                if (validation >= test_result[index]):
                    rank = rankdata(test_result, method="average")
                    rank_list.append(len(gene_set_) - rank[index])
                    break

        validation_rank[disease] = rank_list

    top_n_acc=[top_k(validation_rank,1),top_k(validation_rank,3),top_k(validation_rank,10),top_k(validation_rank,30),top_k(validation_rank,50),top_k(validation_rank,100)]



    average_score=[]
    for disease in validation_rank.keys():
        average_score.extend(validation_rank[disease])
    mean_rank=np.mean(average_score)


    return validation_rank,mean_rank,top_n_acc

def test(model_,evaluation_data,embed_dic_,gene_set_,ranking_score):
    model_.eval()
    validation_rank = dict()
    for disease in evaluation_data:
        positive_genes=disease_genes[disease]

        disease_vec=look_up_embed(embed_dic_,disease)


        testing_list = []
        for gene in gene_set_:
            gene_vec = look_up_embed(embed_dic_,gene)
            testing_list.append([gene_vec,disease_vec])

        validation_list = []
        for gene in positive_genes:
            gene_vec = look_up_embed(embed_dic_,gene)
            validation_list.append([gene_vec,disease_vec])

        testing_list=convert_to_torch_format(testing_list)
        validation_list=convert_to_torch_format(validation_list)
        testing_list=Variable(testing_list)
        validation_list=Variable(validation_list)

        test_result=sorted(model_.predict(testing_list).data,reverse=True)
        valid_result=model_.predict(validation_list).data

        rank_list = []


        for validation in valid_result:

            for index in range(len(test_result)):

                if (validation >= test_result[index]):
                    rank = rankdata(test_result, method="average")
                    rank_list.append(len(gene_set_) - rank[index])
                    break

        validation_rank[disease] = rank_list

    top_n_acc=[top_k(validation_rank,1),top_k(validation_rank,3),top_k(validation_rank,10),top_k(validation_rank,30),top_k(validation_rank,50),top_k(validation_rank,100)]

    average_score=[]


    for disease in validation_rank.keys():
        average_score.extend(validation_rank[disease])
        for gene,ranking_gene in zip(disease_genes[disease],validation_rank[disease]) :
            disease_ge = disease+"_"+gene

            ranking_score[disease_ge]= int(ranking_gene)
    mean_rank=np.mean(average_score)


    return validation_rank,mean_rank,top_n_acc,ranking_score


def count_h(k,rank_list):
    if (k==0):
        return 0
    threshold = k / len(gene_list)
    count = 0
    for value in rank_list:
        if value <= threshold:
            count += 1

    return count / len(rank_list)


def ranked_auc(rank_list):
    rank_dic = {}
    for i in range(1, len(gene_list)):
        rank_dic[i] = count_h(i,rank_list)
    auc = 0
    prior = 10000
    for data in rank_dic.values():
        if (prior == 10000):
            prior = data
        else:
            auc += (1 / 2) * (prior + data) / (len(gene_list) - 1)
            prior = data
    return auc


def calculate_auc(rank_result):
    rank_list = []

    for disease in rank_result.keys():
        for score in rank_result[disease]:
            rank= score/len(gene_list)
            rank_list.append(rank)
    auc=ranked_auc(rank_list)

    return auc




def variance_rank(rank_list):
    highest_rank =0
    lowest_rank =100000
    mean_rank = sum(rank_list)/len(rank_list)
    variance_rank =0

    for data in rank_list:
        if data > highest_rank:
            highest_rank = data
        if data < lowest_rank:
            lowest_rank =data
        variance_rank += (mean_rank- data)*(mean_rank - data)
    variance_rank=variance_rank/len(rank_list)

    return highest_rank, lowest_rank, variance_rank




if __name__ == '__main__':





    with open(associations,"rb") as f:
        disease_genes=pkl.load(f)

    with open(ranking_candidate,"rb") as f:
        gene_list=pkl.load(f)
    gene_list=[gene for gene in gene_list]
    word2vec_model=gensim.models.Word2Vec.load(embedding_model)


    entities=set()
    for disease in disease_genes.keys():
        entities.add(disease)
        for gene in disease_genes[disease]:
            entities.add(gene)
    for gene in gene_list:
        entities.add(gene)

    dic=dict()



    for entity in entities:
        dic[entity]=word2vec_model[entity]




    diseases =[]
    for disease in disease_genes.keys():
        diseases.append(disease)

    random.shuffle(diseases)



    disease_span=int(len(diseases)/10)


    test_datas = [diseases[i:i+disease_span] for i in range(0,len(diseases),disease_span)]

    accumulated_auc = []
    accumulated_rank1=[]
    accumulated_rank3=[]
    accumulated_rank10=[]
    accumulated_rank30=[]
    accumulated_rank50=[]
    accumulated_rank100=[]


    ranking_score=dict()


    for test_data in test_datas:


        g1,d,y,test_disease,train_diseases  =  load_data(dic,disease_genes,gene_list,test_data)
        feature_num=g1.shape[1]

        model=Rank_model(num_feature=feature_num).double()
        opt = optim.Adam(model.parameters(), lr=0.0005, betas=(0.9, 0.999))
        criterion = torch.nn.BCELoss()

        epoches=155

        performance=1000000


        random.shuffle(train_diseases)
        eval_data = train_diseases[int(len(train_diseases)*0.9):]
        training_data = train_diseases[:int(len(train_diseases)*0.9)]

        for epoch in range(epoches):
            g1, d, y= generate_train_data(dic,disease_genes,gene_list,training_data)
            if epoch%10==9:
                auc,rank, top_n_acc=evaluation(model,eval_data,dic,gene_list)

                auc=calculate_auc(auc)
                print("evaluation data set:",auc, rank)
                if rank<performance:
                    torch.save(model, "data/"+str(rank)+"_best_performance.pt")
                    performance=rank


            model,loss_value=train(model,opt,criterion)

        best_model = torch.load("data/"+str(performance)+"_best_performance.pt")

        auc,rank,top_n_acc ,ranking_score= test(best_model, test_data, dic, gene_list,ranking_score)
        auc = calculate_auc(auc)

        print("testing data set:",auc,rank)
        accumulated_auc.append(auc)
        accumulated_rank1.append(top_n_acc[0])
        accumulated_rank3.append(top_n_acc[1])
        accumulated_rank10.append(top_n_acc[2])
        accumulated_rank30.append(top_n_acc[3])
        accumulated_rank50.append(top_n_acc[4])
        accumulated_rank100.append(top_n_acc[5])

    print("the accumulated auc",accumulated_auc)
    highest_auc , lowest_auc, variance_auc = variance_rank(accumulated_auc)
    print("highest auc",highest_auc)
    print("lowest auc", lowest_auc)
    print("difference of auc",highest_auc-lowest_auc)
    print("mean of auc", str(sum(accumulated_auc)/len(accumulated_auc)))
    print("variance of auc", variance_auc)
    print("top 1 accuracy", str(sum(accumulated_rank1)/len(accumulated_rank1)))
    print("top 3 accuracy", str(sum(accumulated_rank3)/len(accumulated_rank3)))
    print("top 5 accuracy", str(sum(accumulated_rank10)/len(accumulated_rank10)))
    print("top 10 accuracy", str(sum(accumulated_rank30)/len(accumulated_rank30)))
    print("top 30 accuracy", str(sum(accumulated_rank50)/len(accumulated_rank50)))
    print("top 100 accuracy", str(sum(accumulated_rank100)/len(accumulated_rank100)))
