# associate the human genes
# associate the mp
# associate the go
# associate the uberon


# define a general method that could generate each type of association


# associate the human disease

# select the gene and disease combinations


# generate the features for the gene and disease

import numpy as np
import pandas as pd
import pickle as pkl


options=["union","intersection","go","uberon","mgi"]

# generate the human disease and hp association
dis_phe=dict()
with open("data/phenotype_annotation.tab","r") as f:
    for line in f.readlines():
        data=line.split("\t")
        if (data[0]=="OMIM")&(data[5][:4]=="OMIM"):
            try:
                dis_phe[data[5].strip()].append(data[4].strip())
            except:
                dis_phe[data[5].strip()]=[data[4].strip()]

print("the number of disease phenotypes ",len(dis_phe))


# obtain the mouse gene and mouse disease association
# obtain the human gene and human disease association


mgi_do=pd.read_table("data/MGI_DO.rpt.txt")
human_gene_disease=[]

for index in mgi_do.index:
    if (mgi_do.loc[index,"Common Organism Name"]=="human"):
        disease=mgi_do.loc[index,"OMIM IDs"]
        sub_dis=disease.split("|")
        human_gene=mgi_do.loc[index,"EntrezGene ID"]
        if disease:
            if human_gene:
                for dis in sub_dis:
                    human_gene_disease.append([int(human_gene),dis.strip()])

mouse_gene_disease=dict()
for index in mgi_do.index:
    mouse_gene=mgi_do.loc[index,"Mouse MGI ID"]
    if str(mouse_gene)!="nan":
        disease=str(mgi_do.loc[index,"OMIM IDs"]).strip()
        sub_dis=disease.split("|")
        for dis in sub_dis:
            if dis !="nan":
                try:
                    mouse_gene_disease[mouse_gene].append(dis.strip())
                except:
                    mouse_gene_disease[mouse_gene]=[dis.strip()]

## generate the dictionary from human to mouse and from mouse to human
human_to_mouse=dict()
mouse_to_human=dict()
geneName_to_id=dict()
with open("data/HMD_HumanPhenotype.rpt.txt",'r') as f:
    for line in f.readlines():
        data=line.split("\t")
        gene_name=data[0].strip()
        human_id=data[1].strip()
        mouse_id=data[5].strip()
        human_to_mouse[human_id]=mouse_id
        mouse_to_human[mouse_id]=human_id
        geneName_to_id[gene_name]=human_id



# obtain the human gene and gene function association


gene_go_feature=dict()
gene_expression_name=set()
with open("data/goa_human.gaf","r") as f:
    for line in f.readlines():
        data=line.split("\t")
        gene_name=data[2].strip()
        evidence_score=data[6].strip()
        go_id=data[4].strip()
        if gene_name in geneName_to_id.keys():
            if not ((evidence_score=="IEA") or (evidence_score=="ND")):
                human_gene=geneName_to_id[gene_name]
                if (human_gene in mouse_to_human.values()):

                    try:
                        gene_go_feature[human_gene].append(go_id)
                    except:
                        gene_go_feature[human_gene]=[go_id]


# obtain the mammalian phenotypes and gene association

mp_pheno=pd.read_table("data/MGI_GenePheno.rpt.txt",names=["Allelic_Composition","Allele_Symbol","Allele_ID","Genetic_Background","Mammalian_Phenotype_ID","PubMed_ID","MGI_Marker","MGI_Genotype_Accession"])
gene_mp_feature=dict()
for index in mp_pheno.index:
    mouse_gene=mp_pheno.loc[index,"MGI_Marker"].strip()
    mouse_pheno=mp_pheno.loc[index,"Mammalian_Phenotype_ID"].strip()
    if (mouse_gene in mouse_to_human.keys()):
        human_gene=str(int(mouse_to_human[mouse_gene]))
        try:
            gene_mp_feature[human_gene].append(mouse_pheno)
        except:
            gene_mp_feature[human_gene]=[mouse_pheno]

# obtain the gene and uberon associations
gene_expression=pd.read_table("data/E-MTAB-5214-query-results.tpms-2.tsv")
gene_expression=gene_expression.drop(columns="EBV-transformed lymphocyte")
gene_expression=gene_expression.drop(columns="transformed skin fibroblast")

anatomy_uberon={'Brodmann (1909) area 24': 'UBERON_0006101',
 'Brodmann (1909) area 9': 'UBERON_0013540',
 'C1 segment of cervical spinal cord': 'UBERON_0006469',
 'adrenal gland': 'UBERON_0002369',
 'amygdala': 'UBERON_0001876',
 'aorta': 'UBERON_0000947',
 'atrium auricular region': 'UBERON_0006618',
 'breast': 'UBERON_0000310',
 'caudate nucleus': 'UBERON_0001873',
 'cerebellar hemisphere': 'UBERON_0002245',
 'cerebral cortex': 'UBERON_0000956',
 'coronary artery': 'UBERON_0001621',
 'cortex of kidney': 'UBERON_0001225',
 'ectocervix': 'UBERON_0012249',
 'endocervix': 'UBERON_0000458',
 'esophagogastric junction': 'UBERON_0007650',
 'esophagus mucosa': 'UBERON_0002469',
 'esophagus muscularis mucosa': 'UBERON_0004648',
 'fallopian tube': 'UBERON_0003889',
 'greater omentum': 'UBERON_0005448',
 'heart left ventricle': 'UBERON_0002084',
 'hypothalamus': 'UBERON_0001898',
 'lower leg skin': 'UBERON_0004264',
 'minor salivary gland': 'UBERON_0001830',
 'nucleus accumbens': 'UBERON_0001882',
 'ovary': 'UBERON_0000992',
 'pancreas': 'UBERON_0001264',
 'pituitary gland': 'UBERON_0000007',
 'prostate gland': 'UBERON_0002367',
 'putamen': 'UBERON_0001874',
 'sigmoid colon': 'UBERON_0001159',
 "small intestine Peyer's patch": 'UBERON_0003454',
 'stomach': 'UBERON_0000945',
 'subcutaneous adipose tissue': 'UBERON_0002190',
 'substantia nigra': 'UBERON_0002038',
 'suprapubic skin': 'UBERON_0036149',
 'testis': 'UBERON_0000473',
 'thyroid gland': 'UBERON_0002046',
 'tibial artery': 'UBERON_0007610',
 'tibial nerve': 'UBERON_0001323',
 'transverse colon': 'UBERON_0001157',
 'urinary bladder': 'UBERON_0001255',
 'uterus': 'UBERON_0000995',
 'vagina': 'UBERON_0000996',
  'blood':"UBERON_0000178",
  'liver':"UBERON_0002107",
  'lung':"UBERON_0002048",
  "spleen":"UBERON_0002106",
  "cerebellum":"UBERON_0002037",
  "skeletal muscle tissue":"UBERON_0001134",
  "hippocampus proper":"UBERON_0002305"}


# convert the nan value into 0 value
for index in gene_expression.index:
    for name in gene_expression.columns[2:]:
        if (str(gene_expression.loc[index,name])=="nan"):
            gene_expression.loc[index,name]=0


# now filter the uberon that has expression value more than 4.0

gene_uberon_feature=dict()
threshold=4.0
gene_uberon_feature_0=dict()
gene_uberon_feature_1_0=dict()
gene_uberon_feature_2_0=dict()
gene_uberon_feature_3_0=dict()
gene_uberon_feature_4_0=dict()
gene_uberon_feature_5_0=dict()
gene_uberon_feature_6_0=dict()
threshold=[0,1.0,2.0,3.0,4.0,5.0,6.0]

for th in threshold:
    for index in gene_expression.index:
        name=gene_expression.loc[index,"Gene Name"]
        if name in geneName_to_id.keys():
            name=geneName_to_id[name]
            if name in human_to_mouse.keys():
                for column in gene_expression.columns[2:]:
                    if gene_expression.loc[index,column]>=th:
                        try:
                            if th==0:
                                gene_uberon_feature_0[name].add(anatomy_uberon[column])
                            elif th==1.0:
                                gene_uberon_feature_1_0[name].add(anatomy_uberon[column])
                            elif th==2.0:
                                gene_uberon_feature_2_0[name].add(anatomy_uberon[column])
                            elif th==3.0:
                                gene_uberon_feature_3_0[name].add(anatomy_uberon[column])
                            elif th==4.0:
                                gene_uberon_feature_4_0[name].add(anatomy_uberon[column])
                            elif th==5.0:
                                gene_uberon_feature_5_0[name].add(anatomy_uberon[column])
                            elif th==6.0:
                                gene_uberon_feature_6_0[name].add(anatomy_uberon[column])
                        except:
                            temp_set=set()
                            temp_set.add(anatomy_uberon[column])
                            # gene_uberon_feature[name]=temp_set
                            if th==0:
                                gene_uberon_feature_0[name]=temp_set
                            elif th==1.0:
                                gene_uberon_feature_1_0[name]=temp_set
                            elif th==2.0:
                                gene_uberon_feature_2_0[name]=temp_set
                            elif th==3.0:
                                gene_uberon_feature_3_0[name]=temp_set
                            elif th==4.0:
                                gene_uberon_feature_4_0[name]=temp_set
                            elif th==5.0:
                                gene_uberon_feature_5_0[name]=temp_set
                            elif th==6.0:
                                gene_uberon_feature_6_0[name]=temp_set

"""
now we have the data set of uberon, mammalian phenotypes and gene function
gene_go_feature
gene_mp_feature
gene_uberon_feature

"""




# generate the intersection data
gene_intersection_feature=dict()

##  updated one
gene_mp_intersection_feature=dict()
gene_go_intersection_feature=dict()
gene_uberon_intersection_feature=dict()
##  updated one


for data in gene_go_feature.keys():
    if data in gene_mp_feature.keys():
        if data in gene_uberon_feature.keys():
            features=set()
            go_intersection_features=set()
            mp_intersection_features=set()
            uberon_intersection_features=set()
            # print(gene_go_feature[data])
            for value in gene_go_feature[data]:
                features.add(value)
            for value in gene_mp_feature[data]:
                features.add(value)
            for value in gene_uberon_feature[data]:
                features.add(value)

# updated features
            for value in gene_go_feature[data]:
                go_intersection_features.add(value)
            for value in gene_mp_feature[data]:
                mp_intersection_features.add(value)
            for value in gene_uberon_feature[data]:
                uberon_intersection_features.add(value)

            gene_go_intersection_feature[data]=go_intersection_features
            gene_mp_intersection_feature[data]=mp_intersection_features
            gene_uberon_intersection_feature[data]=uberon_intersection_features
# updated features

            gene_intersection_feature[data]=features

# generate the union data
gene_total=set()
for data in gene_go_feature.keys():
    gene_total.add(data)
for data in gene_uberon_feature.keys():
    gene_total.add(data)
for data in gene_mp_feature.keys():
    gene_total.add(data)

gene_union_feature=dict()
for data in gene_total:
    features = set()
    if data in gene_go_feature.keys():
        for value in gene_go_feature[data]:
            features.add(value)
    if data in gene_mp_feature.keys():
        for value in gene_mp_feature[data]:
            features.add(value)
    if data in gene_uberon_feature.keys():
        for value in gene_uberon_feature[data]:
            features.add(value)
    gene_union_feature[data] = features

print("num of mp  ",str(len(gene_mp_feature)))
print("num of go  ",str(len(gene_go_feature)))
print("num of uberon   ",str(len(gene_uberon_feature)))
print("num of intersection    ",str(len(gene_intersection_feature)))
print("number of go in intersection",str(len(gene_go_intersection_feature)))
print("number of mp in intersection",str(len(gene_mp_intersection_feature)))
print("number of uberon in intersection",str(len(gene_uberon_intersection_feature)))
print("num of union     ",str(len(gene_union_feature)))


'''
generate the gene and disease associations based on mouse gene and mouse diseases
associate the features with gene and diseases
write the gene and disease out according to the specified data format
'''



def generate_features(mouse_gene_disease_association,gene_feature,disease_feature,type):
    file=open("data/"+type+"_association.txt","w")
    disease_gene=dict()
    count=0

    total_gene=set()
    for gene in gene_feature.keys():
        total_gene.add(gene)
        features=gene_feature[gene]
        for feature in features:
            feature=feature.replace(":","_")
            file.write(gene+" "+'<http://purl.obolibrary.org/obo/'+str(feature)+">"+"\n")


    for disease in disease_feature.keys():
        features=disease_feature[disease]
        for feature in features:
            feature = feature.replace(":", "_")
            file.write(disease+" "+'<http://purl.obolibrary.org/obo/'+str(feature)+">"+"\n")

    for key in mouse_gene_disease_association.keys():
        try:
            human_gene=mouse_to_human[key]
            # g_features=gene_feature[human_gene]
            diseases=mouse_gene_disease_association[key]
            for dis in diseases:
                if dis in disease_feature.keys():
                    if human_gene in gene_feature.keys():
                        count+=1
                        # d_features=disease_feature[dis]
                        # for data in d_features:
                        #     file.write(dis+' '+"<http://purl.obolibrary.org/obo/"+str(data)+">"+"\n")
                        # for data in g_features:
                        #     file.write(human_gene+" "+'<http://purl.obolibrary.org/obo/'+str(data)+">"+"\n")

                        try:
                            disease_gene[dis].append(human_gene)
                        except:
                            disease_gene[dis]=[human_gene]
        except:
            pass

    with open("data/"+type+"_disease_gene.pkl","wb") as f:
        pkl.dump(disease_gene,f)
    with open("data/"+type+"_gene_set.pkl","wb") as f:
        pkl.dump(total_gene,f)
    print("the num of gene disease association :",str(count))
    print(" the number of total gene",len(total_gene))
    file.close()


generate_features(mouse_gene_disease,gene_mp_feature,dis_phe,"mp")
generate_features(mouse_gene_disease,gene_go_feature,dis_phe,"go")
generate_features(mouse_gene_disease,gene_uberon_feature,dis_phe,"uberon")
generate_features(mouse_gene_disease,gene_intersection_feature,dis_phe,"intersection")
generate_features(mouse_gene_disease,gene_go_intersection_feature, dis_phe,"go_intersection")
generate_features(mouse_gene_disease, gene_mp_intersection_feature,dis_phe,"mp_intersection")
generate_features(mouse_gene_disease,gene_uberon_intersection_feature,dis_phe, "uberon_intersection")
generate_features(mouse_gene_disease,gene_union_feature,dis_phe,"union")
