# Apply DL2Vec on disease gene prioritization task
We design DL2Vec and apply this method on disease gene prioritization task. 
Here we put the dataset link and data preprocessing file here.
In order to allow other people to reproduce the same performance, we describe more details in the following.

##  Package dependency
please refer the package dependencies as same as the dependencies of DL2vec

## Data set downloading
create a dictionary and put all the following dataset into data folder<br>
- [Genotypes and Mammalian Phenotype Annotations for Marker Type Genes excluding conditional mutation](http://www.informatics.jax.org/downloads/reports/MGI_GenePheno.rpt)
- [Human/Mouse Orthology with phenotypes](http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt)<br>
[Association of genes with Diseases](http://www.informatics.jax.org/downloads/reports/MGI_DO.rpt)<br>
- [phenotype Associations](http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa)<br>
- [RNA-seq from 53 human tissue samples from the Genotype-Tissue Expression (GTEx)](https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv)<br>
- [gene ontology annotations](http://geneontology.org/gene-associations/goa_human.gaf.gz)<br>
- [phenomeNET](http://aber-owl.net/ontology/PhenomeNET/#/Download)<br>



## Data preprocessing
    python data_preprocessing.py
## Data preprocessing with adding ppi data
Download the preprocessed files that map ontologies to proteins into data folder from: [PPI_data](https://bio2vec.cbrc.kaust.edu.sa/data/DL2vec/)<br>


    python data_preprocessing_adding_ppi.py

## Generate the embedding
    python runDL2vec.py -ontology "ontology file" -associations "association_file" -outfile "embedding output file" -entity_list "entities list need generating embeddding"
  
You can find the generated embedding for all the experiments in [DL2vec_Embeddings](https://bio2vec.cbrc.kaust.edu.sa/data/DL2vec/) 

## Train pointwise ranking model 
    python pointwise_ranking.py -associations "disease and gene association file (dictionary format)" -ranking_candidate "all the genes that need to be prioritized (list format) -embedding_model "DL2vec embedding model)

## Paper citation
If you use our model in any terms of academic use, please cite our paper: Chen, Jun, Azza Th Althagafi, and Robert Hoehndorf. "Predicting candidate genes from phenotypes, functions, and anatomical site of expression." (2020) 
