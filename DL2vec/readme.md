# DL2vec

DL2vec is a project that can convert different types of axioms into graph representation, and then generate an embedding representation for each node and edge.

The main covnertion tool is in DL2vec_embed folder

## pre-installled necessary packages
python >= 3.4 <br>
pandas >=0.24.2 <br>
numpy >= 1.16.2 <br>
torch >= 1.0.1 <br>
networkx >= 2.3 <br>
scipy >= 1.2.1 <br>
scikit-learn >=0.20.3 <br>
Groovy (Groovy Version: 2.4.10 JVM: 1.8.0_121) with Grape for dependency management (http://docs.groovy-lang.org/latest/html/documentation/grape.html).

## running DL2vec
1. download all the files from this repository and put it into a directory named with DL2vec
2. generate the embedding representation for each entity
### generate the embeddings
    python runDL2vec.py -ontology "ontology file" -associations "association_file" -outfile "embedding output file" -entity_list "entities list need generating embeddding [optional]"
    
data format examples:

ontology file: phenomeNET ontology
association_file: 

23370 <http://purl.obolibrary.org/obo/GO_0005886> <br>
54502 <http://purl.obolibrary.org/obo/GO_0005634> <br>
54502 <http://purl.obolibrary.org/obo/GO_0003729> <br>
54502 <http://purl.obolibrary.org/obo/GO_0003723> <br>
79989 <http://purl.obolibrary.org/obo/GO_0120170> <br>
79989 <http://purl.obolibrary.org/obo/GO_0036064> <br>
where first term is the gene id and second term is the GO annotations where can be found in the phenomeNET ontlogy <br>

output: the embedding output location <br>

the followings are mandatory Arguments: <br>
1. ontology_file: ontology file congtains ontology in owl format
2. association file: file contains entity-class associations
3. outfile: output file contains the embedding model
If one of these two mandatory files is missing, an error message will be displayed <br>
you can also specify the following optional arguments:<br>
1. -h, --help show the help message and exit
2. window_size: window size for Word2Vec
3. -mincount minimum count value for Word2Vec
4. entity_list: the entity file in which each entity need to start the random walk and generate the embedding

## Output
The script will save a model that can generate embeddings for each entity <br>

## Reference
If you find our work useful, please cite:
## Acknowledgement
there are two scripts [ProcessOntology.groovy](https://github.com/bio-ontology-research-group/DL2Vec/blob/master/DL2vec/ProcessOntology.groovy) and [getMetadata.groovy](https://github.com/bio-ontology-research-group/DL2Vec/blob/master/DL2vec/getMetadata.groovy) that are adapted from [OPA2vec](https://github.com/bio-ontology-research-group/opa2vec).
