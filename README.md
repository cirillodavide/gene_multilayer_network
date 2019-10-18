# A gene-gene multilayer network

----
## Data sources
**Molecular associations**

Molecular associations are retrieved from BioGRID database. Edge attributes are "physical" and "genetic" interactions.
>BioGRID release 3.5.177 [16/10/2019]
    
    perl bin/BioGRID_graph.pl

**Drug-target associations**

Drug-target interactions are retrieved from KEGG BRITE *Target-based Classification of Compounds* database. Edge attributes are drug identifiers.
>KEGG BRITE release br08010 [16/10/2019]

    python3 bin/KEGG_drugs_graph.py

**Variant-disease associations**

Variant-disease associations are retrieved from the Monarch Disease Ontology (MONDO). Edge attributes are disease identifiers.
>MONDO release 2019-09-30 [16/10/2019]

    python3 bin/MONDO_diseases_graph.py

**Pathways associations**

Pathways associations are retrieved from the lowest level pathway diagram of Reactome. Edge attributes are pathways identifiers.
>Reactome release 70 [16/10/2019]

    perl bin/Reactome_pathways_graph.pl

**Metabolic reaction associations**

Metabolic reaction associations are retrieved from Recon3D via BiGG Models. Edge attributes are metabolites identifiers.
>BiGG last update 2019-09-12 [16/10/2019]

    python3 bin/Recon3D_metabolites_graph.py

----

