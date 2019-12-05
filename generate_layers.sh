#generate the molecular associations layer
perl bin/BioGRID_interactions_graph.pl
#generate the drug-target associations layer
python3 bin/KEGG_drugs_graph.py
#generate the variant-disease associations layer
python3 bin/MONDO_diseases_graph.py
#generate the pathway associations layer
perl bin/Reactome_pathways_graph.pl
#generate the metabolic reaction associations layer
python3 bin/Recon3D_metabolites_graph.py
