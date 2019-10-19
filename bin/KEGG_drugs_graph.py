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
import random
import mygene
import urllib
import datetime

now = datetime.datetime.now()

# check latest version of KEGG target-based classification of compounds (https://www.genome.jp/kegg/brite.html)
# Download json from KEGG target-based classification of compounds and parse genes and targets excluding the 'Unclassified' category.

drug_dict = {}
kegg_dict = {}
release = 'br08310'

def compile_dicts(x0,x1):
	match = re.search(r'\[(.*?)\]',x0['name'])
	if match is not None:
		kegg_gene = match.group(1)
		hsa_match = re.match(r'^HSA:',kegg_gene) # HSA_VAR, for example, are expluded (there is only one)
		if hsa_match is not None:
			kegg_gene = re.sub('HSA:','',kegg_gene).split(' ')
			drug = x1['name']
			drug_id = drug.split(" ")[0]
			drug = re.sub(drug_id,'',drug)
			kegg_dict[drug_id] = kegg_gene
			drug_dict[drug_id] = re.sub(' \(.*?$','',drug)
		else:
			pass
	else:
		pass

url = 'https://www.genome.jp/kegg-bin/download_htext?htext='+release+'.keg&format=json&filedir='
filedata = urllib.request.urlopen(url)
datatowrite = filedata.read()
with open('src/'+release+'.json', 'wb') as f:
	f.write(datatowrite)

with open('src/'+release+'.json') as f:
    d0 = json.load(f)

for d1 in d0['children']:
	for d2 in d1['children']:
		if d2['name'] != 'Unclassified':
			for d3 in d2['children']:
				if 'children' in d3:
					for d4 in d3['children']:
						if 'children' in d4:
							for d5 in d4['children']:
								compile_dicts(d4,d5)

# Connect genes by shared drug (one-to-one gene-target associations will be exluded):

G = nx.Graph()
for k,v in kegg_dict.items():
    G.add_edges_from(itertools.combinations(v, 2), label=k)

nx.write_edgelist(G,"networks/KEGG_drugs."+now.strftime("%d-%m-%Y")+".gr",data=['label'])
