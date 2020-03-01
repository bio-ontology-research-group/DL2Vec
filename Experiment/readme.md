# apply DL2Vec on disease gene prioritization task
We design DL2Vec and apply this method on disease gene prioritization task. 
Here we put the dataset link and data preprocessing file here.
In order to allow other people to reproduce the same performance, we describe more details in the following.

## Data set downloading
create a dictionary and put the all the following dataset into data folder<br>
[Human/Mouse Orthology with phenotypes](http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt)<br>
[Association of genes with Diseases](http://www.informatics.jax.org/downloads/reports/MGI_DO.rpt)<br>
[phenotype Associations](https://hpo.jax.org/app/download/annotation)<br>
[RNA-seq from 53 human tissue samples from the Genotype-Tissue Expression (GTEx)](https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5214/resources/ExperimentDownloadSupplier.RnaSeqBaseline/tpms.tsv)<br>
[gene ontology annotations](http://geneontology.org/gene-associations/goa_human.gaf.gz)<br>
[phenomeNET](http://aber-owl.net/ontology/PhenomeNET/#/Download)<br>


## data preprocessing
  python datapreprocessing.py

## generate the embedding


## train pointwise ranking model 

## Paper citation

