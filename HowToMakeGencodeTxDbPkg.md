# How to make TxDb for genome buildhg38
The aim of this document is to show how to make Bioconductor Transcripts database (TxDb) and package for Gencode track. Generally there are many ways to do it, and here I like to show how to extrack the current Gencode track from UCSC genome browser, which might be tricky, especially when Bionconductor 

## Build Gencode TxDb From UCSC Genome Brower
Let's start with load two major libraries and open a browser session on R workspace.


```r
library(rtracklayer)
library(GenomicFeatures)
session <- browserSession()
genome(session) <- "hg38"
```

To check the track names and table names of Gencode, you can do

```r
rtracklayer::trackNames(session)
```

Currently on on UCSC Genome Brower (9/26/2016), the track name  for Gencode is _Gencdoe v24_ and the table name is _knownGene_.  If you have tried


```r
txdb <- makeTxDbFromUCSC(genome="hg38", tablename="knownGene")
```
and if it works. Then you can ignore the rest.

If it doesn't work, that means there are some errors inhibit in the *GenomicFeatures* pakcage. You will need to take a detour and can use the following code as an example.


```r
tablename <- "knownGene"
track <- "GENCODE v24"
ucsc_txtable <-
    getTable(ucscTableQuery(session, track=track, table=tablename))
mapdef <- GenomicFeatures:::.howToGetTxName2GeneIdMapping(tablename)
txname2geneid <- GenomicFeatures:::.fetchTxName2GeneIdMappingFromUCSC(session, 
           track, tablename, mapdef)
transcript_ids = NULL
txdb <- GenomicFeatures:::.makeTxDbFromUCSCTxTable(ucsc_txtable,
        txname2geneid$genes, 
        genome="hg38", tablename, track, txname2geneid$gene_id_type, 
        full_dataset = is.null(transcript_ids), circ_seqs = DEFAULT_CIRC_SEQS,
        taxonomyId = NA, miRBaseBuild = NA)
```

To build a TxDb package,

```r
makeTxDbPackage(txdb, version = "1.24.0",
                maintainer = "Chao-Jen Wong <cwon2@fredhutch.org>",
                author = "Chao-Jen Wong", destDir = "~/tapscott/hg38")
```
				
## Ensembl From BioMart
Using the most recent annotation on Ensembl? 

```r
library(biomaRt)
library(GenomicFeatures)
listMarts(host="www.ensembl.org")
datasets <- listDatasets(biomaRt::useMart(
    biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org"))
ensembl.txdb <- makeTxDbFromBiomart(biomart="ENSEMBL_MART_ENSEMBL",
                                    dataset="hsapiens_gene_ensembl",
                                    host="www.ensembl.org")
makeTxDbPackage(ensembl.txdb, version="1.85.0", author="Chao-Jen Wong",
                maintainer="Chao-Jen Wong <cwon2@fredhutch.org>")
```
				
## Session Information

```r
sessionInfo()
```

```
## R Under development (unstable) (2016-09-25 r71358)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.3 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.24.5 AnnotationDbi_1.34.4   Biobase_2.32.0        
##  [4] rtracklayer_1.32.2     GenomicRanges_1.24.2   GenomeInfoDb_1.8.7    
##  [7] IRanges_2.6.1          S4Vectors_0.10.3       BiocGenerics_0.18.0   
## [10] knitr_1.14             BiocInstaller_1.22.3  
## 
## loaded via a namespace (and not attached):
##  [1] XVector_0.12.1             magrittr_1.5              
##  [3] zlibbioc_1.18.0            GenomicAlignments_1.8.4   
##  [5] BiocParallel_1.6.6         stringr_1.1.0             
##  [7] tools_3.4.0                SummarizedExperiment_1.2.3
##  [9] DBI_0.5-1                  formatR_1.4               
## [11] bitops_1.0-6               biomaRt_2.28.0            
## [13] RCurl_1.95-4.8             evaluate_0.9              
## [15] RSQLite_1.0.0              stringi_1.1.1             
## [17] Biostrings_2.40.2          Rsamtools_1.24.0          
## [19] XML_3.98-1.4
```
