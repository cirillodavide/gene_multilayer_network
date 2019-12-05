# A gene-gene multilayer network from public resources

### Motivation

To populate the relational graph with molecular information, we created a gene-gene multilayer network consisting of five layers in which nodes represent genes and edges represent different types of associations retrieved from publicly available knowledge bases. The gene-gene multilayer network allows to expand the information of the meduloblastoma data-driven network provided by CURIE upon several layers of molecular information. All the data was downloaded on October 19, 2019. The scripts to generate the gene-gene multilayer network are available at https://github.com/cirillodavide/gene_multilayer_network.

### Data sources

We retrieved information covering five types of gene associations from publicly available databases. The five layers are presented as edgelists in which pairs of genes (Entrez identifiers) are reported with specific relation subtypes.

##### Molecular associations

In this layer, two genes are connected if a physical or genetic association exists. Molecular associations between human genes were obtained from [BioGRID](https://thebiogrid.org), release 3.5.177. BioGRID is a multi-species database of interactions, curated from high-throughput datasets and individual studies [^fn1]. The relation subtypes of the molecular associations layer are labeled as 'physical' or 'genetic'.

##### Drug-target associations

In this layer, two genes are connected if they are both targets of the same drug. Drug-taget associations between human genes were obteined from [KEGG BRITE *Target-based Classification of Compounds*](https://www.genome.jp/kegg-bin/get_htext?br08010.keg), release br08310. KEGG BRITE [^fn2] is a manually curated database of functional hierarchies of various biological objects, such as Drug classifications. The Target-based Classification of Compounds consists of six categories (Protein-coupled receptors, Nuclear receptors, Ion channels, Transportes, Enzymes, Others). One-to-one and unclassified gene-target associations were excluded. The relation subtypes of the drug-target associations are labeled using KEGG drug identifiers (e.g. D09692).

##### Variant-disease associations

In this layer, two genes are connected if they are both reported to be associated to the same disease in genome-wide association studies (GWAS). Variant-disease associations between human genes were obtained from [Monarch Disease Ontology (MonDO)](https://monarchinitiative.org), release 2019-09-30. MonDO [^fn3] is a multi-species ontology generated by merging and harmonizing multiple disease resources (ORDO/Orphanet, DO, OMIM, MESH, etc.). Gene-disease associations are inferred by integrating gene variants (SNPs, SNVs, QTLs, CNVs, among others) from GWAS data above a certain level of significance. We retrieved MonDO entries with associated OMIM identifiers from the [OWL file](http://purl.obolibrary.org/obo/mondo.owl), filtering for evidence code ECO:0000220 (sequencing assay evidence) through the Monarch Solr search service. The relation subtypes of the variant-disease associations are labeled using MonDO disease identifiers (e.g. MONDO_0007179).

##### Pathway associations

In this layer, two genes are connected if they are both annotated to the same pathway. Pathway associations between human genes were obtained from [Reactome](https://reactome.org), release 70. Reactome [^fn4] is a manually curated pathway database. Associations were retrieved from the lowest level pathway diagram of Reactome hierarchy. We found that all annotations are associated to IEA (inferred from electronic annotations) and TAS (traceable author statement) evidence codes. The relation subtypes of the pathway associations are labeled using Reactome pathway identifiers (e.g. R-HSA-72163).

##### Metabolic reaction associations

In this layer, two genes are connected if they are involved in metabolic reactions where product metabolites of one reaction are reactant metabolites of the other one. Metabolic reaction associations between human genes were obtained from Recon3D [^fn5] through [BiGG Models](http://bigg.ucsd.edu), release 2019-09-12. Recon3D is the largest human metabolic network model. Superconnected metabolites (e.g. ATP, CO~2~, H~2~O) [^fn6] were excluded. The relation subtypes of the metabolic reaction associations are labeled using standard metabolite identifiers (e.g. 13dampp).

### Summary

|||||||
|---|---|---|---|---|---|
|relation type|Molecular associations|Drug-target associations|Variant-disease associations|Pathways associations|Metabolic reaction associations|
|relation subtypes|physical, genetic|D09692, D09692,...|MONDO_0007179, MONDO_0000208,...|R-HSA-72163, R-HSA-983712,...|13dampp, 13_cis_retn,...|
|number of nodes|17794|453|976|10718|1786|
|number of edges|343120|1187|16537|888805|52077|

[^fn1]: Oughtred et al. The BioGRID interaction database: 2019 update. Nucleic Acids Res. 2019 Jan 8;47(D1):D529-D541. doi: 10.1093/nar/gky1079. PMID: 30476227
[^fn2]: Kanehisa et al. New approach for understanding genome variations in KEGG. Nucleic Acids Res. 2019 Jan 8;47(D1):D590-D595. doi: 10.1093/nar/gky962. PMID: 30321428
[^fn3]: Mungall et al. The Monarch Initiative: an integrative data and analytic platform connecting phenotypes to genotypes across species. Nucleic Acids Res. 2017 Jan 4;45(D1):D712-D722. doi: 10.1093/nar/gkw1128. PMID: 27899636
[^fn4]: Fabregat et al. The Reactome Pathway Knowledgebase. Nucleic Acids Res. 2018 Jan 4;46(D1):D649-D655. doi: 10.1093/nar/gkx1132. PMID: 29145629
[^fn5]: Brunk et al. Recon3D enables a three-dimensional view of gene variation in human metabolism. Nat Biotechnol. 2018 Mar;36(3):272-281. doi: 10.1038/nbt.4072. PMID: 29457794
[^fn6]: Croes at al. Inferring meaningful pathways in weighted metabolic networks. J Mol Biol. 2006 Feb 10;356(1):222-36. PMID: 16337962
