import rdflib
import re
import pandas as pd
from collections import defaultdict, Counter
from joblib import Parallel, delayed
import multiprocessing
import numpy as np
from pathlib import Path
import json
from pandas.io.json import json_normalize
import networkx as nx
import itertools
from networkx.drawing.nx_agraph import graphviz_layout
import random
import mygene
import urllib
import datetime
import math

now = datetime.datetime.now()

# We use rdflib to extract MONDO ids (and associated OMIM ids) from the latest owl file. We then use the collected MONDO ids to retrieve the associated genes with evidence code ECO:0000220 (sequencing assay evidence). Both steps are time consuming, so we skip them if the output file has already been generated.

file = Path("src/mondo2genes.txt")

if file.exists():
    
    print( 'Variants already retrieved!' )

else:
    
    g=rdflib.Graph()
    g.load('http://purl.obolibrary.org/obo/mondo.owl') # load the ontology

    rdfl = rdflib.Namespace('http://www.w3.org/2000/01/rdf-schema#')
    owl = rdflib.Namespace('http://www.w3.org/2002/07/owl#')
    oboInOwl = rdflib.Namespace('http://www.geneontology.org/formats/oboInOwl#')

    mondo2omim = {}
    mondo2descr = {}
    for subj in g.subjects():
        for i in g.triples((subj,owl['annotatedSource'],None)):
            mondo = i[2]
            mondoID = mondo.toPython().rsplit('/', 1)[-1]
            for j in g.triples((mondo,rdfl['label'],None)):
                descr = j[2].toPython()
                mondo2descr[mondoID] = descr
            for w in g.triples((mondo,oboInOwl['hasDbXref'],None)):
                if re.match('OMIM:',w[2]):
                    omim = w[2].toPython()
                    mondo2omim[mondoID] = omim

    def genes(url):
        c = pd.read_csv(url,sep='\t',header=0).dropna(subset=['evidence'])
        c = c[c['subject_taxon_label'].str.contains('Homo sapiens', na=False)]
        c['evidence'] = c['evidence'].astype('object')
        genes = c[c['evidence'].str.contains('ECO:0000220')].subject_label.tolist()
        genes = list(set(genes))
        return genes

    inputs = []
    for mondo in mondo2omim.keys():
        mondo = mondo.replace('_',':')
        url = 'https://solr.monarchinitiative.org/solr/golr/select/?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=subject,subject_label,subject_taxon,subject_taxon_label,object,object_label,relation,relation_label,evidence,evidence_label,source,is_defined_by,qualifier&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&facet.method=enum&csv.encapsulator=%22&csv.separator=%09&csv.header=true&csv.mv.separator=%7C&fq=subject_category:%22gene%22&fq=object_category:%22disease%22&fq=object_closure:%22'+str(mondo)+'%22&facet.field=subject_taxon_label&q=*:*'
        inputs.append(url)

    num_cores = multiprocessing.cpu_count()
    results = Parallel(n_jobs=num_cores)(delayed(genes)(i) for i in inputs)

    with open("src/mondo2genes.txt","w") as f:
        for (a,b,c) in zip(mondo2omim.keys(),[mondo2descr[x] for x in mondo2omim.keys()],results):
            f.write("{0},{1},{2}\n".format(a,b,c))


# create the graph

df = pd.read_csv('src/mondo2genes.txt', delimiter=',(?=\[)', header=None, names = ['MONDO','genes'], engine='python')
df[['MONDO','disease']] = df['MONDO'].str.split(',',1,expand=True)
df = df[~df['genes'].isin(['[]'])]
df['genes'] = df['genes'].map(lambda x: re.sub('\[|\]|\'','',x))
df['genes'][df['genes'].str.contains('RefSeq')] =  df['genes'][df['genes'].str.contains('RefSeq')].str.replace('\s.*\(RefSeq\)','').str.replace('^\w+ ','')
df['genes'] = df['genes'].map(lambda x: re.sub('\s+|\(human\)', '', x))

df = pd.DataFrame(df.genes.str.split(',').tolist(), index=df.MONDO).stack()
df = df.reset_index()[[0, 'MONDO']]
df.columns = ['genes', 'MONDO']

mg = mygene.MyGeneInfo()
eg = mg.querymany(df['genes'].drop_duplicates().tolist(),scopes='symbol',species=9606,as_dataframe=True)
d = eg.groupby('symbol')['entrezgene'].apply(list).to_dict()
f = eg.groupby('symbol')['_score'].apply(list).to_dict()
for k in f.keys():
	if math.isnan(float(d[k][0])): # if the highest mygene score is nan, use the second highest scoring entrez id (ATXN8OS), else pass
		try:
			d[k] = d[k][1]
		except:
			pass
lst = []
for i in df['genes']:
    try:
        lst.append(d[i][0])
    except:
        lst.append(i) # not found in mygene (C2orf71 and HGNC:12731)
df['entrezgene'] = lst
mondo_dict = df.groupby('MONDO')['entrezgene'].apply(list).to_dict()

G = nx.Graph()
for k,v in mondo_dict.items():
    G.add_edges_from(itertools.combinations(v, 2), label=k) # disease with one gene will be excluded

nx.write_edgelist(G,"networks/MoNDO_diseases."+now.strftime("%d-%m-%Y")+".gr",data=['label'])
