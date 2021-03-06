---
title: Annotation
subtitle: Functional annotation and more 
author: Konstantin Shakhbazov
github: {user: kn3in, repo: annotation, branch: "gh-pages"}
framework: minimal
mode: selfcontained
widgets: [polycharts]
highlighter: highlight.js
hitheme: solarized_light
assets:
  css: 
    - "http://fonts.googleapis.com/css?family=Open+Sans"
    - "http://fonts.googleapis.com/css?family=Open+Sans+Condensed:700"
---



# Intro
We consider various mapping/relationship between gene, SNP, pathway etc and other biological annotations in order to get a biological insight, whatever it may mean in your particular experiment.

- - -

# Genome Browsers and alike
Easiest way to annotate a small list of genes, SNPs, microarray probes etc is to query various web
based interfaces often provided within genome browsers. Also useful when you need quick sanity check of your annotation scripts.

General interfaces
* [UCSC genome browser][]
    * [UCSC table browser][]
    * [ENCODE at UCSC][]
* [Ensembl genome browser](http://www.ensembl.org/Homo_sapiens/Info/Index)
    * [ENCODE at Ensembl][]
* [NCBI](http://www.ncbi.nlm.nih.gov/)
    * [GTEx eQTL browser](http://www.ncbi.nlm.nih.gov/gtex/GTEX2/gtex.cgi)
    * [NCBI Gene P53 example](http://www.ncbi.nlm.nih.gov/gene/7157)
    * [NCBI SNP rs4783244 example][]
* [Biomart][]
* [Gencode][] `human transcriptome annotation; here for the sake of completeness, part of the ENCODE project. Available as a track in the UCSC browser`

[UCSC genome browser]: http://genome.ucsc.edu/cgi-bin/hgGateway
[UCSC table browser]: http://genome.ucsc.edu/cgi-bin/hgTables?command=start
[ENCODE at UCSC]: http://genome.ucsc.edu/ENCODE/
[ENCODE at Ensembl]: http://www.ensembl.org/info/website/tutorials/encode.html
[Gencode]: http://www.gencodegenes.org/
[NCBI SNP rs4783244 example]: http://www.ncbi.nlm.nih.gov/projects/SNP/snp_ref.cgi?rs=4783244
[Gemini]: http://gemini.readthedocs.org/en/latest/index.html
[Browser interface]: http://gemini.readthedocs.org/en/latest/content/browser.html
[Biomart]: http://biomart.org/

SNP
* [Haploreg2](http://www.broadinstitute.org/mammals/haploreg/haploreg.php)
* [Gemini][] `standalone tool`
    * [Browser interface][]
    * [How to use plink files with gemini][]

[How to use plink files with gemini]: http://gemini.readthedocs.org/en/latest/content/faq.html#how-can-i-use-plink-files-with-gemini

Gene Expression
* [BioGPS](http://biogps.org/ "Tissue specific expression")

Pathway Analysis
* [DAVID](http://david.abcc.ncifcrf.gov/ "aka Dreaded pathway analysis")

Biomart, UCSC, Ensembl and NCBI integrate broad spectrum of data within a unified interface, on the other hand there are plenty of annotation tools addressing particular needs e.g. eQTL browsers, miRNA and lincRNA databases not(yet) mentioned here.

- - -
# Biomart
from [www.biomart.org](http://www.biomart.org)
>BioMart is a freely available, open source, federated database system that provides unified access to disparate, geographically distributed
>data sources. It is designed to be data agnostic and platform independent, such that existing databases can easily be incorporated into the
>BioMart framework.

Ensembl is one of the databases available through Biomart. There are numerous APIs to access Ensembl programmatically. We're going to use R/Bioconductor interface which unlike other APIs allows access entire Biomart.

Ensembl/Biomart APIs:
* web interface [Mart view][]
* R/Bioconductor [biomaRt](http://www.bioconductor.org/packages/release/bioc/html/biomaRt.html)
* python: [Pycogent](http://pycogent.org/cookbook/accessing_databases.html#ensembl) `does much more than just querying Ensembl, please look through the documentation`
* ruby: [ruby-ensembl-api](https://github.com/jandot/ruby-ensembl-api)
* perl: [Perl API](http://www.ensembl.org/info/docs/api/index.html)
* language agnostic beta [REST API](http://beta.rest.ensemblgenomes.org/documentation/user_guide)

[Mart view]: http://www.ensembl.org/biomart/martview/

It is instructive to retrieve some
data using the `highest` level web-interface and the `lowest` MySQL level to see that after all
there is no magic. The rest of APIs provide nice in-between level in your preferred language which
hides all particularities of db schema but allows you to deal with the annotation programmatically.

#### Web interface a.k.a. build a query by selecting drop-down menus.
Hopefully this is self explanatory:

0. Go to  [Mart view][]
1. CHOOSE DATABASE: Ensembl Variation 73
2. CHOOSE DATASET: Homo Sapience Short Variation
3. Filters: GENERAL VARIATION FILTERS: Filter by Variation Name: 
 - rs7775397
 - rs4783244
 - rs6450176
4. Attributes: SEQUENCE VARIATION: Variation Name, Minor allele (ALL), 1000 genomes global MAF (ALL)
5. Results

By the same token one can build a query via Biomart:
[result](http://central.biomart.org/martwizard/#!/Genome?mart=Ensembl+72+Short+Variations+\(SNPs+and+indels\)&datasets=hsapiens_snp&step=4&variation_source=dbSNP&snp_filter=rs7775397%2Crs4783244%2Crs6450176&attributes=refsnp_id%2Cchrom_start%2Cchr_name%2Cminor_allele%2Cminor_allele_freq)


#### Direct MySQL query
Query MySQL backend (not recommended, plus no sequence retrieval; schema available [here](http://www.ensembl.org/info/docs/api/index.html))


```bash
mysql --host=ensembldb.ensembl.org -P 3306  --user=anonymous \
-A -e "SELECT v.name, v.minor_allele, v.minor_allele_freq \
FROM variation v WHERE v.name IN ('rs7775397', 'rs4783244', 'rs6450176');" \
homo_sapiens_variation_73_37
```

```
name	minor_allele	minor_allele_freq
rs4783244	T	0.3393
rs6450176	A	0.3356
rs7775397	G	0.033
```

#### Note

Before we go into discovering lots of goodies in the [Bioconductor](http://www.bioconductor.org) project. Bioconductor [conference and tutorial](http://www.bioconductor.org/help/events/) materials
are available [here](http://www.bioconductor.org/help/course-materials/). Each package in the project has a vignette which is worth checking.
Bioconductor also have workflows section describing particular usage of the project e.g. [Using bioconductor for annotation.](http://www.bioconductor.org/help/workflows/annotation/annotation/)

#### biomaRt


```r
library(biomaRt)
```

### Available marts
Contrast ensembl versions available through web/MySQL interfaces and biomaRt.

```r
kable(head(listMarts()))
```

<table>
 <thead>
  <tr>
   <th align="left"> biomart </th>
   <th align="left"> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> ensembl </td>
   <td align="left"> ENSEMBL GENES 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> snp </td>
   <td align="left"> ENSEMBL VARIATION 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> functional_genomics </td>
   <td align="left"> ENSEMBL REGULATION 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> vega </td>
   <td align="left"> VEGA 53  (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> fungi_mart_21 </td>
   <td align="left"> ENSEMBL FUNGI 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> fungi_variations_21 </td>
   <td align="left"> ENSEMBL FUNGI VARIATION 21 (EBI UK) </td>
  </tr>
</tbody>
</table>

### Available datasets

Select `snp` mart and see which datasets are available for the mart.


```r
mart <- useMart("snp")
kable(head(listDatasets(mart)))
```

<table>
 <thead>
  <tr>
   <th align="left"> dataset </th>
   <th align="left"> description </th>
   <th align="left"> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> pabelii_snp </td>
   <td align="left"> Pongo abelii Short Variation (SNPs and indels) (PPYG2) </td>
   <td align="left"> PPYG2 </td>
  </tr>
  <tr>
   <td align="left"> ecaballus_snp </td>
   <td align="left"> Equus caballus Short Variation (SNPs and indels) (EquCab2) </td>
   <td align="left"> EquCab2 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_snp </td>
   <td align="left"> Homo sapiens Short Variation (SNPs and indels) (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_structvar </td>
   <td align="left"> Homo sapiens Structural Variation (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> oanatinus_snp </td>
   <td align="left"> Ornithorhynchus anatinus Short Variation (SNPs and indels) (OANA5) </td>
   <td align="left"> OANA5 </td>
  </tr>
  <tr>
   <td align="left"> tnigroviridis_snp </td>
   <td align="left"> Tetraodon nigroviridis Short Variation (SNPs and indels) (TETRAODON8.0) </td>
   <td align="left"> TETRAODON8.0 </td>
  </tr>
</tbody>
</table>

Select dataset with human snps `hsapiens_snp`.

```r
snpmart <- useDataset("hsapiens_snp", mart = mart)
```

Shortcut in case you know which mart and dataset you are after:

```r
snpmart <- useMart("snp", dataset = "hsapiens_snp")
```
### Attributes and Filters

A dataset queried on a set of fields: `Attributes` i.e. the desired output.
A query narrowed based on a set of fields: `Filters`.

Here are `Attributes` for the hsapiens_snp dataset of the snp mart.

```r
kable(head(listAttributes(snpmart)))
```

<table>
 <thead>
  <tr>
   <th align="left"> name </th>
   <th align="left"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> refsnp_id </td>
   <td align="left"> Variation Name </td>
  </tr>
  <tr>
   <td align="left"> refsnp_source </td>
   <td align="left"> Variation source </td>
  </tr>
  <tr>
   <td align="left"> refsnp_source_description </td>
   <td align="left"> Variation source description </td>
  </tr>
  <tr>
   <td align="left"> chr_name </td>
   <td align="left"> Chromosome name </td>
  </tr>
  <tr>
   <td align="left"> chrom_start </td>
   <td align="left"> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td align="left"> chrom_strand </td>
   <td align="left"> Strand </td>
  </tr>
</tbody>
</table>

Here are the `Filters`

```r
kable(head(listFilters(snpmart)))
```

<table>
 <thead>
  <tr>
   <th align="left"> name </th>
   <th align="left"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> chr_name </td>
   <td align="left"> Chromosome name </td>
  </tr>
  <tr>
   <td align="left"> chrom_start </td>
   <td align="left"> Start </td>
  </tr>
  <tr>
   <td align="left"> chrom_end </td>
   <td align="left"> End </td>
  </tr>
  <tr>
   <td align="left"> band_start </td>
   <td align="left"> Band Start </td>
  </tr>
  <tr>
   <td align="left"> band_end </td>
   <td align="left"> Band End </td>
  </tr>
  <tr>
   <td align="left"> marker_end </td>
   <td align="left"> Marker End </td>
  </tr>
</tbody>
</table>

Example query:
Suppose we have a list of rs SNP ids (`Filter`, we want only data for those rs ids) and would like to figure out where those SNPs located, their minor alleles and MAFs (`Attributes`, we need only those fields returned from biomart). The `getBM` is the main function to query Biomart. We have already seen first four arguments: `attributes`, `filters`, `value` (of the filters) and `mart`.


```r
?getBM
```


```r
top_ids <- c("rs7775397", "rs4783244", "rs6450176")
snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start",
                                 "minor_allele", "minor_allele_freq"),
                    filters = c("snp_filter"),
                      value = top_ids,
                       mart = snpmart)
kable(snp_pos)
```

<table>
 <thead>
  <tr>
   <th align="left"> refsnp_id </th>
   <th align="left"> chr_name </th>
   <th align="right"> chrom_start </th>
   <th align="left"> minor_allele </th>
   <th align="right"> minor_allele_freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> rs4783244 </td>
   <td align="left"> 16 </td>
   <td align="right"> 82662268 </td>
   <td align="left"> T </td>
   <td align="right"> 0.3393 </td>
  </tr>
  <tr>
   <td align="left"> rs6450176 </td>
   <td align="left"> 5 </td>
   <td align="right"> 53298025 </td>
   <td align="left"> A </td>
   <td align="right"> 0.3361 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> 6 </td>
   <td align="right"> 32261252 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> HSCHR6_MHC_MANN </td>
   <td align="right"> 32300691 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> HSCHR6_MHC_COX </td>
   <td align="right"> 32209826 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> HSCHR6_MHC_QBL </td>
   <td align="right"> 32219074 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> HSCHR6_MHC_DBB </td>
   <td align="right"> 32237102 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
  <tr>
   <td align="left"> rs7775397 </td>
   <td align="left"> HSCHR6_MHC_SSTO </td>
   <td align="right"> 32268189 </td>
   <td align="left"> G </td>
   <td align="right"> 0.0326 </td>
  </tr>
</tbody>
</table>

#### All together

```r
library(biomaRt)
snpmart <- useMart("snp", dataset = "hsapiens_snp")
top_ids <- c("rs7775397", "rs4783244", "rs6450176")
snp_pos <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start",
                                 "minor_allele", "minor_allele_freq"),
                    filters = c("snp_filter"),
                      value = top_ids,
                       mart = snpmart)
```

- - -
# Full tables for biomaRt

### Marts

```r
kable(listMarts())
```

<table>
 <thead>
  <tr>
   <th align="left"> biomart </th>
   <th align="left"> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> ensembl </td>
   <td align="left"> ENSEMBL GENES 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> snp </td>
   <td align="left"> ENSEMBL VARIATION 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> functional_genomics </td>
   <td align="left"> ENSEMBL REGULATION 75 (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> vega </td>
   <td align="left"> VEGA 53  (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> fungi_mart_21 </td>
   <td align="left"> ENSEMBL FUNGI 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> fungi_variations_21 </td>
   <td align="left"> ENSEMBL FUNGI VARIATION 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> metazoa_mart_21 </td>
   <td align="left"> ENSEMBL METAZOA 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> metazoa_variations_21 </td>
   <td align="left"> ENSEMBL METAZOA VARIATION 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> plants_mart_21 </td>
   <td align="left"> ENSEMBL PLANTS 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> plants_variations_21 </td>
   <td align="left"> ENSEMBL PLANTS VARIATION 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> protists_mart_21 </td>
   <td align="left"> ENSEMBL PROTISTS 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> protists_variations_21 </td>
   <td align="left"> ENSEMBL PROTISTS VARIATION 21 (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> msd </td>
   <td align="left"> MSD (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> htgt </td>
   <td align="left"> WTSI MOUSE GENETICS PROJECT (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> REACTOME </td>
   <td align="left"> REACTOME (CSHL US) </td>
  </tr>
  <tr>
   <td align="left"> WS220 </td>
   <td align="left"> WORMBASE 220 (CSHL US) </td>
  </tr>
  <tr>
   <td align="left"> biomart </td>
   <td align="left"> MGI (JACKSON LABORATORY US) </td>
  </tr>
  <tr>
   <td align="left"> pride </td>
   <td align="left"> PRIDE (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> prod-intermart_1 </td>
   <td align="left"> INTERPRO (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> unimart </td>
   <td align="left"> UNIPROT (EBI UK) </td>
  </tr>
  <tr>
   <td align="left"> biomartDB </td>
   <td align="left"> PARAMECIUM GENOME (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td align="left"> biblioDB </td>
   <td align="left"> PARAMECIUM BIBLIOGRAPHY (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td align="left"> Eurexpress Biomart </td>
   <td align="left"> EUREXPRESS (MRC EDINBURGH UK) </td>
  </tr>
  <tr>
   <td align="left"> phytozome_mart </td>
   <td align="left"> PHYTOZOME (JGI/CIG US) </td>
  </tr>
  <tr>
   <td align="left"> HapMap_rel27 </td>
   <td align="left"> HAPMAP 27 (NCBI US) </td>
  </tr>
  <tr>
   <td align="left"> CosmicMart </td>
   <td align="left"> COSMIC (SANGER UK) </td>
  </tr>
  <tr>
   <td align="left"> cildb_all_v2 </td>
   <td align="left"> CILDB INPARANOID AND FILTERED BEST HIT (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td align="left"> cildb_inp_v2 </td>
   <td align="left"> CILDB INPARANOID (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td align="left"> experiments </td>
   <td align="left"> INTOGEN EXPERIMENTS </td>
  </tr>
  <tr>
   <td align="left"> oncomodules </td>
   <td align="left"> INTOGEN ONCOMODULES </td>
  </tr>
  <tr>
   <td align="left"> gmap_japonica </td>
   <td align="left"> RICE-MAP JAPONICA (PEKING UNIVESITY CHINA) </td>
  </tr>
  <tr>
   <td align="left"> europhenomeannotations </td>
   <td align="left"> EUROPHENOME </td>
  </tr>
  <tr>
   <td align="left"> ikmc </td>
   <td align="left"> IKMC GENES AND PRODUCTS (IKMC) </td>
  </tr>
  <tr>
   <td align="left"> EMAGE gene expression </td>
   <td align="left"> EMAGE GENE EXPRESSION </td>
  </tr>
  <tr>
   <td align="left"> EMAP anatomy ontology </td>
   <td align="left"> EMAP ANATOMY ONTOLOGY </td>
  </tr>
  <tr>
   <td align="left"> EMAGE browse repository </td>
   <td align="left"> EMAGE BROWSE REPOSITORY </td>
  </tr>
  <tr>
   <td align="left"> GermOnline </td>
   <td align="left"> GERMONLINE </td>
  </tr>
  <tr>
   <td align="left"> Sigenae_Oligo_Annotation_Ensembl_61 </td>
   <td align="left"> SIGENAE OLIGO ANNOTATION (ENSEMBL 61) </td>
  </tr>
  <tr>
   <td align="left"> Sigenae Oligo Annotation (Ensembl 59) </td>
   <td align="left"> SIGENAE OLIGO ANNOTATION (ENSEMBL 59) </td>
  </tr>
  <tr>
   <td align="left"> Sigenae Oligo Annotation (Ensembl 56) </td>
   <td align="left"> SIGENAE OLIGO ANNOTATION (ENSEMBL 56) </td>
  </tr>
  <tr>
   <td align="left"> Breast_mart_69 </td>
   <td align="left"> BCCTB Bioinformatics Portal (UK and Ireland) </td>
  </tr>
  <tr>
   <td align="left"> K562_Gm12878 </td>
   <td align="left"> Predictive models of gene regulation from processed high-throughput epigenomics data: K562 vs. Gm12878 </td>
  </tr>
  <tr>
   <td align="left"> Hsmm_Hmec </td>
   <td align="left"> Predictive models of gene regulation from processed high-throughput epigenomics data: Hsmm vs. Hmec </td>
  </tr>
  <tr>
   <td align="left"> Pancreas63 </td>
   <td align="left"> PANCREATIC EXPRESSION DATABASE (BARTS CANCER INSTITUTE UK) </td>
  </tr>
  <tr>
   <td align="left"> Public_OBIOMART </td>
   <td align="left"> Genetic maps (markers, Qtls), Polymorphisms (snps, genes), Genetic and Phenotype resources with Genes annotations </td>
  </tr>
  <tr>
   <td align="left"> Public_VITIS </td>
   <td align="left"> Grapevine 8x, stuctural annotation with Genetic maps (genetic markers..) </td>
  </tr>
  <tr>
   <td align="left"> Public_VITIS_12x </td>
   <td align="left"> Grapevine 12x, stuctural and functional annotation with Genetic maps (genetic markers..) </td>
  </tr>
  <tr>
   <td align="left"> Prod_WHEAT </td>
   <td align="left"> Wheat, stuctural annotation with Genetic maps (genetic markers..) and Polymorphisms (snps) </td>
  </tr>
  <tr>
   <td align="left"> Public_TAIRV10 </td>
   <td align="left"> Arabidopsis Thaliana TAIRV10, genes functional annotation </td>
  </tr>
  <tr>
   <td align="left"> Public_MAIZE </td>
   <td align="left"> Zea mays ZmB73, genes functional annotation </td>
  </tr>
  <tr>
   <td align="left"> Prod_POPLAR </td>
   <td align="left"> Populus trichocarpa, genes functional annotation </td>
  </tr>
  <tr>
   <td align="left"> Prod_POPLAR_V2 </td>
   <td align="left"> Populus trichocarpa, genes functional annotation V2.0 </td>
  </tr>
  <tr>
   <td align="left"> Prod_BOTRYTISEDIT </td>
   <td align="left"> Botrytis cinerea T4, genes functional annotation  </td>
  </tr>
  <tr>
   <td align="left"> Prod_ </td>
   <td align="left"> Botrytis cinerea B0510, genes functional annotation  </td>
  </tr>
  <tr>
   <td align="left"> Prod_SCLEROEDIT </td>
   <td align="left"> Sclerotinia sclerotiorum, genes functional annotation  </td>
  </tr>
  <tr>
   <td align="left"> Prod_LMACULANSEDIT </td>
   <td align="left"> Leptosphaeria maculans, genes functional annotation </td>
  </tr>
  <tr>
   <td align="left"> vb_mart_22 </td>
   <td align="left"> VectorBase Genes </td>
  </tr>
  <tr>
   <td align="left"> vb_snp_mart_22 </td>
   <td align="left"> VectorBase Variation </td>
  </tr>
  <tr>
   <td align="left"> expression </td>
   <td align="left"> VectorBase Expression </td>
  </tr>
  <tr>
   <td align="left"> ENSEMBL_MART_PLANT </td>
   <td align="left"> GRAMENE 40 ENSEMBL GENES (CSHL/CORNELL US) </td>
  </tr>
  <tr>
   <td align="left"> ENSEMBL_MART_PLANT_SNP </td>
   <td align="left"> GRAMENE 40 VARIATION (CSHL/CORNELL US) </td>
  </tr>
</tbody>
</table>
### Datasets available for the snp mart

```r
kable(listDatasets(mart))
```

<table>
 <thead>
  <tr>
   <th align="left"> dataset </th>
   <th align="left"> description </th>
   <th align="left"> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> pabelii_snp </td>
   <td align="left"> Pongo abelii Short Variation (SNPs and indels) (PPYG2) </td>
   <td align="left"> PPYG2 </td>
  </tr>
  <tr>
   <td align="left"> ecaballus_snp </td>
   <td align="left"> Equus caballus Short Variation (SNPs and indels) (EquCab2) </td>
   <td align="left"> EquCab2 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_snp </td>
   <td align="left"> Homo sapiens Short Variation (SNPs and indels) (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_structvar </td>
   <td align="left"> Homo sapiens Structural Variation (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> oanatinus_snp </td>
   <td align="left"> Ornithorhynchus anatinus Short Variation (SNPs and indels) (OANA5) </td>
   <td align="left"> OANA5 </td>
  </tr>
  <tr>
   <td align="left"> tnigroviridis_snp </td>
   <td align="left"> Tetraodon nigroviridis Short Variation (SNPs and indels) (TETRAODON8.0) </td>
   <td align="left"> TETRAODON8.0 </td>
  </tr>
  <tr>
   <td align="left"> ggallus_snp </td>
   <td align="left"> Gallus gallus Short Variation (SNPs and indels) (Galgal4) </td>
   <td align="left"> Galgal4 </td>
  </tr>
  <tr>
   <td align="left"> oaries_snp </td>
   <td align="left"> Ovis Aries Short Variation (SNPs and indels) (Oar_v3.1) </td>
   <td align="left"> Oar_v3.1 </td>
  </tr>
  <tr>
   <td align="left"> scerevisiae_snp </td>
   <td align="left"> Saccharomyces cerevisiae Short Variation (SNPs and indels) (R64-1-1) </td>
   <td align="left"> R64-1-1 </td>
  </tr>
  <tr>
   <td align="left"> drerio_structvar </td>
   <td align="left"> Danio rerio Structural Variation (Zv9) </td>
   <td align="left"> Zv9 </td>
  </tr>
  <tr>
   <td align="left"> mmulatta_structvar </td>
   <td align="left"> Macaca mulatta Structural Variation (MMUL_1) </td>
   <td align="left"> MMUL_1 </td>
  </tr>
  <tr>
   <td align="left"> mmusculus_snp </td>
   <td align="left"> Mus musculus Short Variation (SNPs and indels) (GRCm38.p2) </td>
   <td align="left"> GRCm38.p2 </td>
  </tr>
  <tr>
   <td align="left"> mmusculus_structvar </td>
   <td align="left"> Mus musculus Structural Variation (GRCm38.p2) </td>
   <td align="left"> GRCm38.p2 </td>
  </tr>
  <tr>
   <td align="left"> drerio_snp </td>
   <td align="left"> Danio rerio Short Variation (SNPs and indels) (Zv9) </td>
   <td align="left"> Zv9 </td>
  </tr>
  <tr>
   <td align="left"> mdomestica_snp </td>
   <td align="left"> Monodelphis domestica Short Variation (SNPs and indels) (monDom5) </td>
   <td align="left"> monDom5 </td>
  </tr>
  <tr>
   <td align="left"> cfamiliaris_structvar </td>
   <td align="left"> Canis familiaris Structural Variation (CanFam3.1) </td>
   <td align="left"> CanFam3.1 </td>
  </tr>
  <tr>
   <td align="left"> btaurus_structvar </td>
   <td align="left"> Bos taurus Structural Variation (UMD3.1) </td>
   <td align="left"> UMD3.1 </td>
  </tr>
  <tr>
   <td align="left"> ptroglodytes_snp </td>
   <td align="left"> Pan troglodytes Short Variation (SNPs and indels) (CHIMP2.1.4) </td>
   <td align="left"> CHIMP2.1.4 </td>
  </tr>
  <tr>
   <td align="left"> btaurus_snp </td>
   <td align="left"> Bos taurus Short Variation (SNPs and indels) (UMD3.1) </td>
   <td align="left"> UMD3.1 </td>
  </tr>
  <tr>
   <td align="left"> mmulatta_snp </td>
   <td align="left"> Macaca mulatta Short Variation (SNPs and indels) (MMUL_1) </td>
   <td align="left"> MMUL_1 </td>
  </tr>
  <tr>
   <td align="left"> nleucogenys_snp </td>
   <td align="left"> Nomascus leucogenys Short Variation (SNPs and indels) (Nleu1.0) </td>
   <td align="left"> Nleu1.0 </td>
  </tr>
  <tr>
   <td align="left"> mgallopavo_snp </td>
   <td align="left"> Meleagris gallopavo Short Variation (SNPs and indels) (UMD2) </td>
   <td align="left"> UMD2 </td>
  </tr>
  <tr>
   <td align="left"> sscrofa_structvar </td>
   <td align="left"> Sus scrofa Structural Variation (Sscrofa10.2) </td>
   <td align="left"> Sscrofa10.2 </td>
  </tr>
  <tr>
   <td align="left"> ecaballus_structvar </td>
   <td align="left"> Equus caballus Structural Variation (EquCab2) </td>
   <td align="left"> EquCab2 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_snp_som </td>
   <td align="left"> Homo sapiens Somatic Short Variation (SNPs and indels) (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> tguttata_snp </td>
   <td align="left"> Taeniopygia guttata Short Variation (SNPs and indels) (taeGut3.2.4) </td>
   <td align="left"> taeGut3.2.4 </td>
  </tr>
  <tr>
   <td align="left"> fcatus_snp </td>
   <td align="left"> Felis catus Short Variation (SNPs and indels) (Felis_catus_6.2) </td>
   <td align="left"> Felis_catus_6.2 </td>
  </tr>
  <tr>
   <td align="left"> cfamiliaris_snp </td>
   <td align="left"> Canis familiaris Short Variation (SNPs and indels) (CanFam3.1) </td>
   <td align="left"> CanFam3.1 </td>
  </tr>
  <tr>
   <td align="left"> hsapiens_structvar_som </td>
   <td align="left"> Homo sapiens Somatic Structural Variation (GRCh37.p13) </td>
   <td align="left"> GRCh37.p13 </td>
  </tr>
  <tr>
   <td align="left"> sscrofa_snp </td>
   <td align="left"> Sus scrofa Short Variation (SNPs and indels) (Sscrofa10.2) </td>
   <td align="left"> Sscrofa10.2 </td>
  </tr>
  <tr>
   <td align="left"> dmelanogaster_snp </td>
   <td align="left"> Drosophila melanogaster Short Variation (SNPs and indels) (BDGP5) </td>
   <td align="left"> BDGP5 </td>
  </tr>
  <tr>
   <td align="left"> rnorvegicus_snp </td>
   <td align="left"> Rattus norvegicus Short Variation (SNPs and indels) (Rnor_5.0) </td>
   <td align="left"> Rnor_5.0 </td>
  </tr>
</tbody>
</table>
### Attributes for the snp mart, the hsapiens_snp dataset

```r
kable(listAttributes(snpmart))
```

<table>
 <thead>
  <tr>
   <th align="left"> name </th>
   <th align="left"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> refsnp_id </td>
   <td align="left"> Variation Name </td>
  </tr>
  <tr>
   <td align="left"> refsnp_source </td>
   <td align="left"> Variation source </td>
  </tr>
  <tr>
   <td align="left"> refsnp_source_description </td>
   <td align="left"> Variation source description </td>
  </tr>
  <tr>
   <td align="left"> chr_name </td>
   <td align="left"> Chromosome name </td>
  </tr>
  <tr>
   <td align="left"> chrom_start </td>
   <td align="left"> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td align="left"> chrom_strand </td>
   <td align="left"> Strand </td>
  </tr>
  <tr>
   <td align="left"> allele </td>
   <td align="left"> Variant Alleles </td>
  </tr>
  <tr>
   <td align="left"> mapweight </td>
   <td align="left"> Mapweight </td>
  </tr>
  <tr>
   <td align="left"> validated </td>
   <td align="left"> Evidence status </td>
  </tr>
  <tr>
   <td align="left"> allele_1 </td>
   <td align="left"> Ancestral allele </td>
  </tr>
  <tr>
   <td align="left"> minor_allele </td>
   <td align="left"> Minor allele (ALL) </td>
  </tr>
  <tr>
   <td align="left"> minor_allele_freq </td>
   <td align="left"> 1000 genomes global MAF (ALL) </td>
  </tr>
  <tr>
   <td align="left"> minor_allele_count </td>
   <td align="left"> 1000 genomes global MAC (ALL) </td>
  </tr>
  <tr>
   <td align="left"> clinical_significance </td>
   <td align="left"> Clinical significance </td>
  </tr>
  <tr>
   <td align="left"> synonym_name </td>
   <td align="left"> Synonym name </td>
  </tr>
  <tr>
   <td align="left"> synonym_source </td>
   <td align="left"> Synonym source </td>
  </tr>
  <tr>
   <td align="left"> synonym_source_description </td>
   <td align="left"> Synonym source description </td>
  </tr>
  <tr>
   <td align="left"> variation_names </td>
   <td align="left"> Associated variation names </td>
  </tr>
  <tr>
   <td align="left"> study_type </td>
   <td align="left"> Study type </td>
  </tr>
  <tr>
   <td align="left"> study_external_ref </td>
   <td align="left"> Study External Reference </td>
  </tr>
  <tr>
   <td align="left"> study_description </td>
   <td align="left"> Study Description </td>
  </tr>
  <tr>
   <td align="left"> source_name </td>
   <td align="left"> Source name </td>
  </tr>
  <tr>
   <td align="left"> associated_gene </td>
   <td align="left"> Associated gene with phenotype </td>
  </tr>
  <tr>
   <td align="left"> phenotype_description </td>
   <td align="left"> Phenotype description </td>
  </tr>
  <tr>
   <td align="left"> phenotype_significance </td>
   <td align="left"> Phenotype significance [0 non significant, 1 significant] </td>
  </tr>
  <tr>
   <td align="left"> associated_variant_risk_allele </td>
   <td align="left"> Associated variant risk allele </td>
  </tr>
  <tr>
   <td align="left"> p_value </td>
   <td align="left"> P value </td>
  </tr>
  <tr>
   <td align="left"> set_name </td>
   <td align="left"> Variation Set Name </td>
  </tr>
  <tr>
   <td align="left"> set_description </td>
   <td align="left"> Variation Set Description </td>
  </tr>
  <tr>
   <td align="left"> title_20137 </td>
   <td align="left"> Title </td>
  </tr>
  <tr>
   <td align="left"> authors_20137 </td>
   <td align="left"> Authors </td>
  </tr>
  <tr>
   <td align="left"> year_20137 </td>
   <td align="left"> Year </td>
  </tr>
  <tr>
   <td align="left"> pmid_20137 </td>
   <td align="left"> PubMed ID </td>
  </tr>
  <tr>
   <td align="left"> pmcid_20137 </td>
   <td align="left"> PMC reference number (PMCID) </td>
  </tr>
  <tr>
   <td align="left"> ucsc_id_20137 </td>
   <td align="left"> UCSC ID </td>
  </tr>
  <tr>
   <td align="left"> doi_20137 </td>
   <td align="left"> Digital Object Identifier </td>
  </tr>
  <tr>
   <td align="left"> ensembl_gene_stable_id </td>
   <td align="left"> Ensembl Gene ID </td>
  </tr>
  <tr>
   <td align="left"> ensembl_transcript_stable_id </td>
   <td align="left"> Ensembl Transcript ID </td>
  </tr>
  <tr>
   <td align="left"> ensembl_transcript_chrom_strand </td>
   <td align="left"> Transcript strand </td>
  </tr>
  <tr>
   <td align="left"> ensembl_type </td>
   <td align="left"> Biotype </td>
  </tr>
  <tr>
   <td align="left"> consequence_type_tv </td>
   <td align="left"> Consequence to transcript </td>
  </tr>
  <tr>
   <td align="left"> consequence_allele_string </td>
   <td align="left"> Consequence specific allele </td>
  </tr>
  <tr>
   <td align="left"> ensembl_peptide_allele </td>
   <td align="left"> Protein allele </td>
  </tr>
  <tr>
   <td align="left"> cdna_start </td>
   <td align="left"> Variation start in cDNA (bp) </td>
  </tr>
  <tr>
   <td align="left"> cdna_end </td>
   <td align="left"> Variation end in cDNA (bp) </td>
  </tr>
  <tr>
   <td align="left"> translation_start </td>
   <td align="left"> Variation start in translation (aa) </td>
  </tr>
  <tr>
   <td align="left"> translation_end </td>
   <td align="left"> Variation end in translation (aa) </td>
  </tr>
  <tr>
   <td align="left"> cds_start </td>
   <td align="left"> Variation start in CDS (bp) </td>
  </tr>
  <tr>
   <td align="left"> cds_end </td>
   <td align="left"> Variation end in CDS (bp) </td>
  </tr>
  <tr>
   <td align="left"> distance_to_transcript </td>
   <td align="left"> Distance to transcript </td>
  </tr>
  <tr>
   <td align="left"> polyphen_prediction </td>
   <td align="left"> PolyPhen prediction </td>
  </tr>
  <tr>
   <td align="left"> polyphen_score </td>
   <td align="left"> PolyPhen score </td>
  </tr>
  <tr>
   <td align="left"> sift_prediction </td>
   <td align="left"> SIFT prediction </td>
  </tr>
  <tr>
   <td align="left"> sift_score </td>
   <td align="left"> SIFT score </td>
  </tr>
  <tr>
   <td align="left"> feature_stable_id_20126 </td>
   <td align="left"> Regulatory Feature Stable ID </td>
  </tr>
  <tr>
   <td align="left"> allele_string_20126 </td>
   <td align="left"> Regulatory Feature Allele String </td>
  </tr>
  <tr>
   <td align="left"> consequence_types_20126 </td>
   <td align="left"> Regulatory Feature Consequence Type </td>
  </tr>
  <tr>
   <td align="left"> feature_stable_id_20125 </td>
   <td align="left"> Motif Feature Stable ID </td>
  </tr>
  <tr>
   <td align="left"> allele_string_20125 </td>
   <td align="left"> Motif Feature Allele String </td>
  </tr>
  <tr>
   <td align="left"> consequence_types_20125 </td>
   <td align="left"> Motif Feature Consequence Type </td>
  </tr>
  <tr>
   <td align="left"> in_informative_position_20125 </td>
   <td align="left"> High Information Position </td>
  </tr>
  <tr>
   <td align="left"> motif_score_delta_20125 </td>
   <td align="left"> Motif Score Change </td>
  </tr>
  <tr>
   <td align="left"> motif_name_20125 </td>
   <td align="left"> Motif Name </td>
  </tr>
  <tr>
   <td align="left"> motif_start_20125 </td>
   <td align="left"> Motif Position </td>
  </tr>
  <tr>
   <td align="left"> snp </td>
   <td align="left"> Variation sequence </td>
  </tr>
  <tr>
   <td align="left"> upstream_flank </td>
   <td align="left"> upstream_flank </td>
  </tr>
  <tr>
   <td align="left"> downstream_flank </td>
   <td align="left"> downstream_flank </td>
  </tr>
  <tr>
   <td align="left"> chr_name </td>
   <td align="left"> Chromosome name </td>
  </tr>
  <tr>
   <td align="left"> chrom_start </td>
   <td align="left"> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td align="left"> chrom_strand </td>
   <td align="left"> Strand </td>
  </tr>
  <tr>
   <td align="left"> refsnp_id </td>
   <td align="left"> Variation Name </td>
  </tr>
  <tr>
   <td align="left"> refsnp_source </td>
   <td align="left"> Variation source </td>
  </tr>
  <tr>
   <td align="left"> allele </td>
   <td align="left"> Variant Alleles </td>
  </tr>
  <tr>
   <td align="left"> validated </td>
   <td align="left"> Evidence status </td>
  </tr>
  <tr>
   <td align="left"> mapweight </td>
   <td align="left"> Mapweight </td>
  </tr>
  <tr>
   <td align="left"> ensembl_peptide_allele </td>
   <td align="left"> Protein allele </td>
  </tr>
</tbody>
</table>
### Filters for the snp mart, the hsapiens_snp dataset

```r
kable(listFilters(snpmart))
```

<table>
 <thead>
  <tr>
   <th align="left"> name </th>
   <th align="left"> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td align="left"> chr_name </td>
   <td align="left"> Chromosome name </td>
  </tr>
  <tr>
   <td align="left"> chrom_start </td>
   <td align="left"> Start </td>
  </tr>
  <tr>
   <td align="left"> chrom_end </td>
   <td align="left"> End </td>
  </tr>
  <tr>
   <td align="left"> band_start </td>
   <td align="left"> Band Start </td>
  </tr>
  <tr>
   <td align="left"> band_end </td>
   <td align="left"> Band End </td>
  </tr>
  <tr>
   <td align="left"> marker_end </td>
   <td align="left"> Marker End </td>
  </tr>
  <tr>
   <td align="left"> marker_start </td>
   <td align="left"> Marker Start </td>
  </tr>
  <tr>
   <td align="left"> chromosomal_region </td>
   <td align="left"> Chromosome Regions (e.g 1:100:10000:-1,1:100000:200000:1) </td>
  </tr>
  <tr>
   <td align="left"> strand </td>
   <td align="left"> Strand </td>
  </tr>
  <tr>
   <td align="left"> variation_source </td>
   <td align="left"> Variation source </td>
  </tr>
  <tr>
   <td align="left"> snp_filter </td>
   <td align="left"> Filter by Variation Name (e.g. rs123, CM000001) </td>
  </tr>
  <tr>
   <td align="left"> variation_synonym_source </td>
   <td align="left"> Variation Synonym source </td>
  </tr>
  <tr>
   <td align="left"> study_type </td>
   <td align="left"> Study type </td>
  </tr>
  <tr>
   <td align="left"> phenotype_description </td>
   <td align="left"> Phenotype description </td>
  </tr>
  <tr>
   <td align="left"> phenotype_significance </td>
   <td align="left"> Phenotype significance </td>
  </tr>
  <tr>
   <td align="left"> variation_set_name </td>
   <td align="left"> Variation Set Name </td>
  </tr>
  <tr>
   <td align="left"> sift_prediction </td>
   <td align="left"> SIFT Prediction </td>
  </tr>
  <tr>
   <td align="left"> sift_score </td>
   <td align="left"> SIFT score <= </td>
  </tr>
  <tr>
   <td align="left"> polyphen_prediction </td>
   <td align="left"> PolyPhen Prediction </td>
  </tr>
  <tr>
   <td align="left"> polyphen_score </td>
   <td align="left"> PolyPhen score >= </td>
  </tr>
  <tr>
   <td align="left"> minor_allele_freq </td>
   <td align="left"> Global minor allele frequency <= </td>
  </tr>
  <tr>
   <td align="left"> minor_allele_freq_second </td>
   <td align="left"> Global minor allele frequency >= </td>
  </tr>
  <tr>
   <td align="left"> clinical_significance </td>
   <td align="left"> Clinical_significance </td>
  </tr>
  <tr>
   <td align="left"> with_validated </td>
   <td align="left"> Variations that have been validated </td>
  </tr>
  <tr>
   <td align="left"> with_variation_citation </td>
   <td align="left"> Variations with citations </td>
  </tr>
  <tr>
   <td align="left"> distance_to_transcript </td>
   <td align="left"> Distance to transcript <= </td>
  </tr>
  <tr>
   <td align="left"> ensembl_gene </td>
   <td align="left"> Ensembl Gene ID(s) [Max 500] </td>
  </tr>
  <tr>
   <td align="left"> so_parent_name </td>
   <td align="left"> Parent term name </td>
  </tr>
  <tr>
   <td align="left"> feature_stable_id </td>
   <td align="left"> Filter by Regulatory Stable ID(s) (e.g. ENSR00001529861) [Max 500 ADVISED] </td>
  </tr>
  <tr>
   <td align="left"> motif_name </td>
   <td align="left"> Motif Name </td>
  </tr>
</tbody>
</table>









