---
title: Annotation
subtitle: Functional annotation and more 
author: Konstantin Shakhbazov
github: {user: kn3in, repo: nothereyet, branch: "gh-pages"}
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
   <th> biomart </th>
   <th> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> ensembl </td>
   <td> ENSEMBL GENES 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> snp </td>
   <td> ENSEMBL VARIATION 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> functional_genomics </td>
   <td> ENSEMBL REGULATION 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> vega </td>
   <td> VEGA 52  (SANGER UK) </td>
  </tr>
  <tr>
   <td> fungi_mart_19 </td>
   <td> ENSEMBL FUNGI 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> fungi_variations_18 </td>
   <td> ENSEMBL FUNGI VARIATION 19 (EBI UK) </td>
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
   <th> dataset </th>
   <th> description </th>
   <th> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> pabelii_snp </td>
   <td> Pongo abelii Short Variation (SNPs and indels) (PPYG2) </td>
   <td> PPYG2 </td>
  </tr>
  <tr>
   <td> ecaballus_snp </td>
   <td> Equus caballus Short Variation (SNPs and indels) (EquCab2) </td>
   <td> EquCab2 </td>
  </tr>
  <tr>
   <td> hsapiens_snp </td>
   <td> Homo sapiens Short Variation (SNPs and indels) (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> hsapiens_structvar </td>
   <td> Homo sapiens Structural Variation (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> oanatinus_snp </td>
   <td> Ornithorhynchus anatinus Short Variation (SNPs and indels) (OANA5) </td>
   <td> OANA5 </td>
  </tr>
  <tr>
   <td> tnigroviridis_snp </td>
   <td> Tetraodon nigroviridis Short Variation (SNPs and indels) (TETRAODON8.0) </td>
   <td> TETRAODON8.0 </td>
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
   <th> name </th>
   <th> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> refsnp_id </td>
   <td> Variation Name </td>
  </tr>
  <tr>
   <td> refsnp_source </td>
   <td> Variation source </td>
  </tr>
  <tr>
   <td> refsnp_source_description </td>
   <td> Variation source description </td>
  </tr>
  <tr>
   <td> chr_name </td>
   <td> Chromosome name </td>
  </tr>
  <tr>
   <td> chrom_start </td>
   <td> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td> chrom_strand </td>
   <td> Strand </td>
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
   <th> name </th>
   <th> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> chr_name </td>
   <td> Chromosome name </td>
  </tr>
  <tr>
   <td> chrom_start </td>
   <td> Start </td>
  </tr>
  <tr>
   <td> chrom_end </td>
   <td> End </td>
  </tr>
  <tr>
   <td> band_start </td>
   <td> Band Start </td>
  </tr>
  <tr>
   <td> band_end </td>
   <td> Band End </td>
  </tr>
  <tr>
   <td> marker_end </td>
   <td> Marker End </td>
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
   <th> refsnp_id </th>
   <th> chr_name </th>
   <th> chrom_start </th>
   <th> minor_allele </th>
   <th> minor_allele_freq </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> rs4783244 </td>
   <td> 16 </td>
   <td> 82662268 </td>
   <td> T </td>
   <td> 0.3393 </td>
  </tr>
  <tr>
   <td> rs6450176 </td>
   <td>  5 </td>
   <td> 53298025 </td>
   <td> A </td>
   <td> 0.3356 </td>
  </tr>
  <tr>
   <td> rs7775397 </td>
   <td>  6 </td>
   <td> 32261252 </td>
   <td> G </td>
   <td> 0.0330 </td>
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
   <th> biomart </th>
   <th> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> ensembl </td>
   <td> ENSEMBL GENES 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> snp </td>
   <td> ENSEMBL VARIATION 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> functional_genomics </td>
   <td> ENSEMBL REGULATION 72 (SANGER UK) </td>
  </tr>
  <tr>
   <td> vega </td>
   <td> VEGA 52  (SANGER UK) </td>
  </tr>
  <tr>
   <td> fungi_mart_19 </td>
   <td> ENSEMBL FUNGI 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> fungi_variations_18 </td>
   <td> ENSEMBL FUNGI VARIATION 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> metazoa_mart_19 </td>
   <td> ENSEMBL METAZOA 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> metazoa_variations_18 </td>
   <td> ENSEMBL METAZOA VARIATION 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> plants_mart_19 </td>
   <td> ENSEMBL PLANTS 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> plants_variations_18 </td>
   <td> ENSEMBL PLANTS VARIATION 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> protists_mart_19 </td>
   <td> ENSEMBL PROTISTS 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> protists_variations_18 </td>
   <td> ENSEMBL PROTISTS VARIATION 19 (EBI UK) </td>
  </tr>
  <tr>
   <td> msd </td>
   <td> MSD (EBI UK) </td>
  </tr>
  <tr>
   <td> htgt </td>
   <td> WTSI MOUSE GENETICS PROJECT (SANGER UK) </td>
  </tr>
  <tr>
   <td> REACTOME </td>
   <td> REACTOME (CSHL US) </td>
  </tr>
  <tr>
   <td> WS220 </td>
   <td> WORMBASE 220 (CSHL US) </td>
  </tr>
  <tr>
   <td> biomart </td>
   <td> MGI (JACKSON LABORATORY US) </td>
  </tr>
  <tr>
   <td> g4public </td>
   <td> HGNC (EBI UK) </td>
  </tr>
  <tr>
   <td> pride </td>
   <td> PRIDE (EBI UK) </td>
  </tr>
  <tr>
   <td> prod-intermart_1 </td>
   <td> INTERPRO (EBI UK) </td>
  </tr>
  <tr>
   <td> unimart </td>
   <td> UNIPROT (EBI UK) </td>
  </tr>
  <tr>
   <td> biomartDB </td>
   <td> PARAMECIUM GENOME (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td> biblioDB </td>
   <td> PARAMECIUM BIBLIOGRAPHY (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td> Eurexpress Biomart </td>
   <td> EUREXPRESS (MRC EDINBURGH UK) </td>
  </tr>
  <tr>
   <td> phytozome_mart </td>
   <td> PHYTOZOME (JGI/CIG US) </td>
  </tr>
  <tr>
   <td> HapMap_rel27 </td>
   <td> HAPMAP 27 (NCBI US) </td>
  </tr>
  <tr>
   <td> CosmicMart </td>
   <td> COSMIC (SANGER UK) </td>
  </tr>
  <tr>
   <td> cildb_all_v2 </td>
   <td> CILDB INPARANOID AND FILTERED BEST HIT (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td> cildb_inp_v2 </td>
   <td> CILDB INPARANOID (CNRS FRANCE) </td>
  </tr>
  <tr>
   <td> experiments </td>
   <td> INTOGEN EXPERIMENTS </td>
  </tr>
  <tr>
   <td> combinations </td>
   <td> INTOGEN COMBINATIONS </td>
  </tr>
  <tr>
   <td> oncomodules </td>
   <td> INTOGEN ONCOMODULES </td>
  </tr>
  <tr>
   <td> gmap_japonica </td>
   <td> RICE-MAP JAPONICA (PEKING UNIVESITY CHINA) </td>
  </tr>
  <tr>
   <td> europhenomeannotations </td>
   <td> EUROPHENOME </td>
  </tr>
  <tr>
   <td> emma_biomart </td>
   <td> THE EUROPEAN MOUSE MUTANT ARCHIVE (EMMA) </td>
  </tr>
  <tr>
   <td> ikmc </td>
   <td> IKMC GENES AND PRODUCTS (IKMC) </td>
  </tr>
  <tr>
   <td> EMAGE gene expression </td>
   <td> EMAGE GENE EXPRESSION </td>
  </tr>
  <tr>
   <td> EMAP anatomy ontology </td>
   <td> EMAP ANATOMY ONTOLOGY </td>
  </tr>
  <tr>
   <td> EMAGE browse repository </td>
   <td> EMAGE BROWSE REPOSITORY </td>
  </tr>
  <tr>
   <td> GermOnline </td>
   <td> GERMONLINE </td>
  </tr>
  <tr>
   <td> Sigenae_Oligo_Annotation_Ensembl_61 </td>
   <td> SIGENAE OLIGO ANNOTATION (ENSEMBL 61) </td>
  </tr>
  <tr>
   <td> Sigenae Oligo Annotation (Ensembl 59) </td>
   <td> SIGENAE OLIGO ANNOTATION (ENSEMBL 59) </td>
  </tr>
  <tr>
   <td> Sigenae Oligo Annotation (Ensembl 56) </td>
   <td> SIGENAE OLIGO ANNOTATION (ENSEMBL 56) </td>
  </tr>
  <tr>
   <td> Breast_mart_58 </td>
   <td> BREAST CANCER CAMPAIGN TISSUE BANK EXPRESSION DATABASE </td>
  </tr>
  <tr>
   <td> K562_Gm12878 </td>
   <td> Predictive models of gene regulation from processed high-throughput epigenomics data: K562 vs. Gm12878 </td>
  </tr>
  <tr>
   <td> Hsmm_Hmec </td>
   <td> Predictive models of gene regulation from processed high-throughput epigenomics data: Hsmm vs. Hmec </td>
  </tr>
  <tr>
   <td> GC_mart </td>
   <td> GWASmart </td>
  </tr>
  <tr>
   <td> UTRMart </td>
   <td> AURA </td>
  </tr>
  <tr>
   <td> Pancreas63 </td>
   <td> PANCREATIC EXPRESSION DATABASE (BARTS CANCER INSTITUTE UK) </td>
  </tr>
  <tr>
   <td> Public_OBIOMART </td>
   <td> Genetic maps (markers, Qtls), Polymorphisms (snps, genes), Genetic and Phenotype resources with Genes annotations </td>
  </tr>
  <tr>
   <td> Public_VITIS </td>
   <td> Grapevine 8x, stuctural annotation with Genetic maps (genetic markers..) </td>
  </tr>
  <tr>
   <td> Public_VITIS_12x </td>
   <td> Grapevine 12x, stuctural and functional annotation with Genetic maps (genetic markers..) </td>
  </tr>
  <tr>
   <td> Prod_WHEAT </td>
   <td> Wheat, stuctural annotation with Genetic maps (genetic markers..) and Polymorphisms (snps) </td>
  </tr>
  <tr>
   <td> Public_TAIRV10 </td>
   <td> Arabidopsis Thaliana TAIRV10, genes functional annotation </td>
  </tr>
  <tr>
   <td> Public_MAIZE </td>
   <td> Zea mays ZmB73, genes functional annotation </td>
  </tr>
  <tr>
   <td> Prod_POPLAR </td>
   <td> Populus trichocarpa, genes functional annotation </td>
  </tr>
  <tr>
   <td> Prod_POPLAR_V2 </td>
   <td> Populus trichocarpa, genes functional annotation V2.0 </td>
  </tr>
  <tr>
   <td> Prod_BOTRYTISEDIT </td>
   <td> Botrytis cinerea T4, genes functional annotation  </td>
  </tr>
  <tr>
   <td> Prod_ </td>
   <td> Botrytis cinerea B0510, genes functional annotation  </td>
  </tr>
  <tr>
   <td> Prod_SCLEROEDIT </td>
   <td> Sclerotinia sclerotiorum, genes functional annotation  </td>
  </tr>
  <tr>
   <td> Prod_LMACULANSEDIT </td>
   <td> Leptosphaeria maculans, genes functional annotation </td>
  </tr>
  <tr>
   <td> GRAMENE_MAP_37 </td>
   <td> GRAMENE 37 MAPPINGS (CSHL/CORNELL US) </td>
  </tr>
  <tr>
   <td> QTL_MART </td>
   <td> GRAMENE 37 QTL DB (CSHL/CORNELL US) </td>
  </tr>
  <tr>
   <td> vb_mart_19 </td>
   <td> Vectorbase Genes </td>
  </tr>
  <tr>
   <td> vb_snp_mart_19 </td>
   <td> Vectorbase Variation </td>
  </tr>
  <tr>
   <td> expression </td>
   <td> Vectorbase Expression Mart </td>
  </tr>
  <tr>
   <td> ENSEMBL_MART_PLANT </td>
   <td> GRAMENE 37 ENSEMBL GENES (CSHL/CORNELL US) </td>
  </tr>
  <tr>
   <td> ENSEMBL_MART_PLANT_SNP </td>
   <td> GRAMENE 37 VARIATION (CSHL/CORNELL US) </td>
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
   <th> dataset </th>
   <th> description </th>
   <th> version </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> pabelii_snp </td>
   <td> Pongo abelii Short Variation (SNPs and indels) (PPYG2) </td>
   <td> PPYG2 </td>
  </tr>
  <tr>
   <td> ecaballus_snp </td>
   <td> Equus caballus Short Variation (SNPs and indels) (EquCab2) </td>
   <td> EquCab2 </td>
  </tr>
  <tr>
   <td> hsapiens_snp </td>
   <td> Homo sapiens Short Variation (SNPs and indels) (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> hsapiens_structvar </td>
   <td> Homo sapiens Structural Variation (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> oanatinus_snp </td>
   <td> Ornithorhynchus anatinus Short Variation (SNPs and indels) (OANA5) </td>
   <td> OANA5 </td>
  </tr>
  <tr>
   <td> tnigroviridis_snp </td>
   <td> Tetraodon nigroviridis Short Variation (SNPs and indels) (TETRAODON8.0) </td>
   <td> TETRAODON8.0 </td>
  </tr>
  <tr>
   <td> ggallus_snp </td>
   <td> Gallus gallus Short Variation (SNPs and indels) (Galgal4) </td>
   <td> Galgal4 </td>
  </tr>
  <tr>
   <td> scerevisiae_snp </td>
   <td> Saccharomyces cerevisiae Short Variation (SNPs and indels) (EF2) </td>
   <td> EF2 </td>
  </tr>
  <tr>
   <td> drerio_structvar </td>
   <td> Danio rerio Structural Variation (Zv9) </td>
   <td> Zv9 </td>
  </tr>
  <tr>
   <td> mmulatta_structvar </td>
   <td> Macaca mulatta Structural Variation (MMUL1.0) </td>
   <td> MMUL1.0 </td>
  </tr>
  <tr>
   <td> mmusculus_snp </td>
   <td> Mus musculus Short Variation (SNPs and indels) (GRCm38.p1) </td>
   <td> GRCm38.p1 </td>
  </tr>
  <tr>
   <td> mmusculus_structvar </td>
   <td> Mus musculus Structural Variation (GRCm38.p1) </td>
   <td> GRCm38.p1 </td>
  </tr>
  <tr>
   <td> drerio_snp </td>
   <td> Danio rerio Short Variation (SNPs and indels) (Zv9) </td>
   <td> Zv9 </td>
  </tr>
  <tr>
   <td> mdomestica_snp </td>
   <td> Monodelphis domestica Short Variation (SNPs and indels) (monDom5) </td>
   <td> monDom5 </td>
  </tr>
  <tr>
   <td> cfamiliaris_structvar </td>
   <td> Canis familiaris Structural Variation (CanFam3.1) </td>
   <td> CanFam3.1 </td>
  </tr>
  <tr>
   <td> ptroglodytes_snp </td>
   <td> Pan troglodytes Short Variation (SNPs and indels) (CHIMP2.1.4) </td>
   <td> CHIMP2.1.4 </td>
  </tr>
  <tr>
   <td> btaurus_structvar </td>
   <td> Bos taurus Structural Variation (UMD3.1) </td>
   <td> UMD3.1 </td>
  </tr>
  <tr>
   <td> btaurus_snp </td>
   <td> Bos taurus Short Variation (SNPs and indels) (UMD3.1) </td>
   <td> UMD3.1 </td>
  </tr>
  <tr>
   <td> mmulatta_snp </td>
   <td> Macaca mulatta Short Variation (SNPs and indels) (MMUL1.0) </td>
   <td> MMUL1.0 </td>
  </tr>
  <tr>
   <td> nleucogenys_snp </td>
   <td> Nomascus leucogenys Short Variation (SNPs and indels) (Nleu1.0) </td>
   <td> Nleu1.0 </td>
  </tr>
  <tr>
   <td> ecaballus_structvar </td>
   <td> Equus caballus Structural Variation (EquCab2) </td>
   <td> EquCab2 </td>
  </tr>
  <tr>
   <td> sscrofa_structvar </td>
   <td> Sus scrofa Structural Variation (Sscrofa10.2) </td>
   <td> Sscrofa10.2 </td>
  </tr>
  <tr>
   <td> hsapiens_snp_som </td>
   <td> Homo sapiens Somatic Short Variation (SNPs and indels) (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> tguttata_snp </td>
   <td> Taeniopygia guttata Short Variation (SNPs and indels) (taeGut3.2.4) </td>
   <td> taeGut3.2.4 </td>
  </tr>
  <tr>
   <td> fcatus_snp </td>
   <td> Felis catus Short Variation (SNPs and indels) (Felis_catus_6.2) </td>
   <td> Felis_catus_6.2 </td>
  </tr>
  <tr>
   <td> cfamiliaris_snp </td>
   <td> Canis familiaris Short Variation (SNPs and indels) (CanFam3.1) </td>
   <td> CanFam3.1 </td>
  </tr>
  <tr>
   <td> hsapiens_structvar_som </td>
   <td> Homo sapiens Somatic Structural Variation (GRCh37.p11) </td>
   <td> GRCh37.p11 </td>
  </tr>
  <tr>
   <td> sscrofa_snp </td>
   <td> Sus scrofa Short Variation (SNPs and indels) (Sscrofa10.2) </td>
   <td> Sscrofa10.2 </td>
  </tr>
  <tr>
   <td> dmelanogaster_snp </td>
   <td> Drosophila melanogaster Short Variation (SNPs and indels) (BDGP5) </td>
   <td> BDGP5 </td>
  </tr>
  <tr>
   <td> rnorvegicus_snp </td>
   <td> Rattus norvegicus Short Variation (SNPs and indels) (Rnor_5.0) </td>
   <td> Rnor_5.0 </td>
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
   <th> name </th>
   <th> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> refsnp_id </td>
   <td> Variation Name </td>
  </tr>
  <tr>
   <td> refsnp_source </td>
   <td> Variation source </td>
  </tr>
  <tr>
   <td> refsnp_source_description </td>
   <td> Variation source description </td>
  </tr>
  <tr>
   <td> chr_name </td>
   <td> Chromosome name </td>
  </tr>
  <tr>
   <td> chrom_start </td>
   <td> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td> chrom_strand </td>
   <td> Strand </td>
  </tr>
  <tr>
   <td> allele </td>
   <td> Variant Alleles </td>
  </tr>
  <tr>
   <td> mapweight </td>
   <td> Mapweight </td>
  </tr>
  <tr>
   <td> validated </td>
   <td> Evidence status </td>
  </tr>
  <tr>
   <td> allele_1 </td>
   <td> Ancestral allele </td>
  </tr>
  <tr>
   <td> minor_allele </td>
   <td> Minor allele (ALL) </td>
  </tr>
  <tr>
   <td> minor_allele_freq </td>
   <td> 1000 genomes global MAF (ALL) </td>
  </tr>
  <tr>
   <td> minor_allele_count </td>
   <td> 1000 genomes global MAC (ALL) </td>
  </tr>
  <tr>
   <td> clinical_significance </td>
   <td> Clinical significance </td>
  </tr>
  <tr>
   <td> synonym_name </td>
   <td> Synonym name </td>
  </tr>
  <tr>
   <td> synonym_source </td>
   <td> Synonym source </td>
  </tr>
  <tr>
   <td> synonym_source_description </td>
   <td> Synonym source description </td>
  </tr>
  <tr>
   <td> variation_names </td>
   <td> Associated variation names </td>
  </tr>
  <tr>
   <td> study_type </td>
   <td> Study type </td>
  </tr>
  <tr>
   <td> study_external_ref </td>
   <td> Study External Reference </td>
  </tr>
  <tr>
   <td> study_description </td>
   <td> Study Description </td>
  </tr>
  <tr>
   <td> source_name </td>
   <td> Source name </td>
  </tr>
  <tr>
   <td> associated_gene </td>
   <td> Associated gene with phenotype </td>
  </tr>
  <tr>
   <td> phenotype_name </td>
   <td> Phenotype name </td>
  </tr>
  <tr>
   <td> phenotype_description </td>
   <td> Phenotype description </td>
  </tr>
  <tr>
   <td> phenotype_significance </td>
   <td> Phenotype significance [0 non significant, 1 significant] </td>
  </tr>
  <tr>
   <td> associated_variant_risk_allele </td>
   <td> Associated variant risk allele </td>
  </tr>
  <tr>
   <td> p_value </td>
   <td> P value </td>
  </tr>
  <tr>
   <td> set_name </td>
   <td> Variation Set Name </td>
  </tr>
  <tr>
   <td> set_description </td>
   <td> Variation Set Description </td>
  </tr>
  <tr>
   <td> title_20137 </td>
   <td> Title </td>
  </tr>
  <tr>
   <td> authors_20137 </td>
   <td> Authors </td>
  </tr>
  <tr>
   <td> pmid_20137 </td>
   <td> PubMed ID </td>
  </tr>
  <tr>
   <td> pmcid_20137 </td>
   <td> PMC reference number (PMCID) </td>
  </tr>
  <tr>
   <td> ensembl_gene_stable_id </td>
   <td> Ensembl Gene ID </td>
  </tr>
  <tr>
   <td> ensembl_transcript_stable_id </td>
   <td> Ensembl Transcript ID </td>
  </tr>
  <tr>
   <td> ensembl_transcript_chrom_strand </td>
   <td> Transcript strand </td>
  </tr>
  <tr>
   <td> ensembl_type </td>
   <td> Biotype </td>
  </tr>
  <tr>
   <td> consequence_type_tv </td>
   <td> Consequence to transcript </td>
  </tr>
  <tr>
   <td> consequence_allele_string </td>
   <td> Consequence specific allele </td>
  </tr>
  <tr>
   <td> ensembl_peptide_allele </td>
   <td> Protein allele </td>
  </tr>
  <tr>
   <td> cdna_start </td>
   <td> Variation start in cDNA (bp) </td>
  </tr>
  <tr>
   <td> cdna_end </td>
   <td> Variation end in cDNA (bp) </td>
  </tr>
  <tr>
   <td> translation_start </td>
   <td> Variation start in translation (aa) </td>
  </tr>
  <tr>
   <td> translation_end </td>
   <td> Variation end in translation (aa) </td>
  </tr>
  <tr>
   <td> cds_start </td>
   <td> Variation start in CDS (bp) </td>
  </tr>
  <tr>
   <td> cds_end </td>
   <td> Variation end in CDS (bp) </td>
  </tr>
  <tr>
   <td> distance_to_transcript </td>
   <td> Distance to transcript </td>
  </tr>
  <tr>
   <td> polyphen_prediction </td>
   <td> PolyPhen prediction </td>
  </tr>
  <tr>
   <td> polyphen_score </td>
   <td> PolyPhen score </td>
  </tr>
  <tr>
   <td> sift_prediction </td>
   <td> SIFT prediction </td>
  </tr>
  <tr>
   <td> sift_score </td>
   <td> SIFT score </td>
  </tr>
  <tr>
   <td> feature_stable_id_20126 </td>
   <td> Regulatory Feature Stable ID </td>
  </tr>
  <tr>
   <td> allele_string_20126 </td>
   <td> Regulatory Feature Allele String </td>
  </tr>
  <tr>
   <td> consequence_types_20126 </td>
   <td> Regulatory Feature Consequence Type </td>
  </tr>
  <tr>
   <td> feature_stable_id_20125 </td>
   <td> Motif Feature Stable ID </td>
  </tr>
  <tr>
   <td> allele_string_20125 </td>
   <td> Motif Feature Allele String </td>
  </tr>
  <tr>
   <td> consequence_types_20125 </td>
   <td> Motif Feature Consequence Type </td>
  </tr>
  <tr>
   <td> in_informative_position_20125 </td>
   <td> High Information Position </td>
  </tr>
  <tr>
   <td> motif_score_delta_20125 </td>
   <td> Motif Score Change </td>
  </tr>
  <tr>
   <td> motif_name_20125 </td>
   <td> Motif Name </td>
  </tr>
  <tr>
   <td> motif_start_20125 </td>
   <td> Motif Position </td>
  </tr>
  <tr>
   <td> snp </td>
   <td> Variation sequence </td>
  </tr>
  <tr>
   <td> upstream_flank </td>
   <td> upstream_flank </td>
  </tr>
  <tr>
   <td> downstream_flank </td>
   <td> downstream_flank </td>
  </tr>
  <tr>
   <td> chr_name </td>
   <td> Chromosome name </td>
  </tr>
  <tr>
   <td> chrom_start </td>
   <td> Position on Chromosome (bp) </td>
  </tr>
  <tr>
   <td> chrom_strand </td>
   <td> Strand </td>
  </tr>
  <tr>
   <td> refsnp_id </td>
   <td> Variation Name </td>
  </tr>
  <tr>
   <td> refsnp_source </td>
   <td> Variation source </td>
  </tr>
  <tr>
   <td> allele </td>
   <td> Variant Alleles </td>
  </tr>
  <tr>
   <td> validated </td>
   <td> Evidence status </td>
  </tr>
  <tr>
   <td> mapweight </td>
   <td> Mapweight </td>
  </tr>
  <tr>
   <td> ensembl_peptide_allele </td>
   <td> Protein allele </td>
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
   <th> name </th>
   <th> description </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td> chr_name </td>
   <td> Chromosome name </td>
  </tr>
  <tr>
   <td> chrom_start </td>
   <td> Start </td>
  </tr>
  <tr>
   <td> chrom_end </td>
   <td> End </td>
  </tr>
  <tr>
   <td> band_start </td>
   <td> Band Start </td>
  </tr>
  <tr>
   <td> band_end </td>
   <td> Band End </td>
  </tr>
  <tr>
   <td> marker_end </td>
   <td> Marker End </td>
  </tr>
  <tr>
   <td> marker_start </td>
   <td> Marker Start </td>
  </tr>
  <tr>
   <td> chromosomal_region </td>
   <td> Chromosome Regions (e.g 1:100:10000:-1,1:100000:200000:1) </td>
  </tr>
  <tr>
   <td> strand </td>
   <td> Strand </td>
  </tr>
  <tr>
   <td> variation_source </td>
   <td> Variation source </td>
  </tr>
  <tr>
   <td> snp_filter </td>
   <td> Filter by Variation Name (e.g. rs123, CM000001) </td>
  </tr>
  <tr>
   <td> variation_synonym_source </td>
   <td> Variation Synonym source </td>
  </tr>
  <tr>
   <td> study_type </td>
   <td> Study type </td>
  </tr>
  <tr>
   <td> phenotype_description </td>
   <td> Phenotype description </td>
  </tr>
  <tr>
   <td> phenotype_significance </td>
   <td> Phenotype significance </td>
  </tr>
  <tr>
   <td> variation_set_name </td>
   <td> Variation Set Name </td>
  </tr>
  <tr>
   <td> sift_prediction </td>
   <td> SIFT Prediction </td>
  </tr>
  <tr>
   <td> sift_score </td>
   <td> SIFT score <= </td>
  </tr>
  <tr>
   <td> polyphen_prediction </td>
   <td> PolyPhen Prediction </td>
  </tr>
  <tr>
   <td> polyphen_score </td>
   <td> PolyPhen score >= </td>
  </tr>
  <tr>
   <td> minor_allele_freq </td>
   <td> Global minor allele frequency <= </td>
  </tr>
  <tr>
   <td> minor_allele_freq_second </td>
   <td> Global minor allele frequency >= </td>
  </tr>
  <tr>
   <td> clinical_significance </td>
   <td> Clinical_significance </td>
  </tr>
  <tr>
   <td> with_validated </td>
   <td> Variations that have been validated </td>
  </tr>
  <tr>
   <td> with_variation_citation </td>
   <td> Variations with citations </td>
  </tr>
  <tr>
   <td> distance_to_transcript </td>
   <td> Distance to transcript <= </td>
  </tr>
  <tr>
   <td> ensembl_gene </td>
   <td> Ensembl Gene ID(s) [Max 500] </td>
  </tr>
  <tr>
   <td> so_parent_name </td>
   <td> Parent term name </td>
  </tr>
  <tr>
   <td> feature_stable_id </td>
   <td> Filter by Regulatory Stable ID(s) (e.g. ENSR00001529861) [Max 500 ADVISED] </td>
  </tr>
  <tr>
   <td> motif_name </td>
   <td> Motif Name </td>
  </tr>
</tbody>
</table>










