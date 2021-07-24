Quiz 4
================
Harlan Gillespie
23/07/2021

## Dependencies

``` r
library(goseq)
library(Biobase)
library(limma)
library(DESeq2)
library(org.Mm.eg.db)
```

## Question 1

When performing gene set analysis it is critical to use the same
annotation as was used in pre-processing steps. Read the paper behind
the Bottomly data set on the ReCount database:
<http://www.ncbi.nlm.nih.gov/pubmed?term=21455293>

Using the paper and the function: supportedGenomes() in the goseq
package can you figure out which of the Mouse genome builds they aligned
the reads to.

## Answer 1

The paper states that “All reads were realigned to the NCBI m37 version
of the mouse genome assembly”

``` r
sp = supportedGenomes()
sp[sp$species =="Mouse",]
```

    ##     db species      date          name
    ## 78 mm9   Mouse Jul. 2007 NCBI Build 37
    ## 79 mm8   Mouse Feb. 2006 NCBI Build 36
    ## 80 mm7   Mouse Aug. 2005 NCBI Build 35
    ##                                                                                               AvailableGeneIDs
    ## 78 acembly,ccdsGene,ensGene,exoniphy,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,xenoRefGene
    ## 79          ccdsGene,ensGene,geneSymbol,geneid,genscan,knownGene,nscanGene,refGene,sgpGene,sibGene,xenoRefGene
    ## 80                                     ensGene,geneSymbol,geneid,genscan,knownGene,refGene,sgpGene,xenoRefGene

This shows 3 possible mouse genomes, mm9 being the only genome which
matches the NCBI m37 description.

## Question 2

Load the Bottomly data with the following code and perform a
differential expression analysis using limma with only the strain
variable as an outcome. How many genes are differentially expressed at
the 5% FDR level using Benjamini-Hochberg correction? What is the gene
identifier of the first gene differentially expressed at this level
(just in order, not the smallest FDR) ? (hint: the featureNames function
may be useful)

``` r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bottomly_eset.RData")
load(file=con)
close(con)
bot = bottomly.eset
pdata_bot=pData(bot)
fdata_bot = featureData(bot)
edata = exprs(bot)
fdata_bot = fdata_bot[rowMeans(edata) > 5]
edata = edata[rowMeans(edata) > 5, ]
edata = log2(edata+1)
```

## Answer 2

first we fit a linear regression model to each gene’s expression data
using the lm.fit function:

``` r
mod = model.matrix(~ pdata_bot$strain)
fit = lmFit(edata, mod)
fit = eBayes(fit)
```

Now we can apply Benjamini-Hochberg FDR correction:

``` r
fit_BH = topTable(fit, number=dim(edata)[1], adjust.method = "BH", p.value = 0.05, sort.by="none")
```

    ## Removing intercept from test coefficients

Now we can find how many genes are differentially expressed at the 5%
FDR level using Benjamini-Hochberg correction and the gene identifier of
the first gene differentially expressed at this level:

``` r
dim(fit_BH)
```

    ## [1] 223   6

``` r
head(fit_BH, 1)
```

    ##                        logFC  AveExpr         t      P.Value   adj.P.Val
    ## ENSMUSG00000000402 -1.222062 4.292471 -4.509076 5.312399e-05 0.004394846
    ##                           B
    ## ENSMUSG00000000402 1.583059

## Question 3

Use the nullp and goseq functions in the goseq package to perform a gene
ontology analysis. What is the top category that comes up as over
represented? (hint: you will need to use the genome information on the
genome from question 1 and the differential expression analysis from
question 2.)

## Answer 3

First, we create an integer list showing the DE (differentially
expressed at the 5% FDR level) status of each gene in the data set:

``` r
genes = as.integer(rownames(edata) %in% rownames(fit_BH))
names(genes) = rownames(edata)
```

Next, we set up a weighting function for the genes using the nullp
function. This is necessary due to the size of genes introducing some
bias into the expression levels so this should be accounted for.

``` r
pwf=nullp(genes,"mm9","ensGene")
```

    ## Loading mm9 length data...

![](Quiz-4_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Here we use a parametric test to look for differences in enrichment with
respect to different categories:

``` r
GO.wall = goseq(pwf, "mm9", "ensGene")
```

    ## Fetching GO annotations...

    ## For 502 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
head(GO.wall)
```

    ##         category over_represented_pvalue under_represented_pvalue numDEInCat
    ## 1945  GO:0004888            3.557394e-07                0.9999999         25
    ## 8998  GO:0038023            7.811515e-07                0.9999998         27
    ## 12558 GO:0060089            7.811515e-07                0.9999998         27
    ## 464   GO:0001653            7.501345e-06                0.9999990         10
    ## 3106  GO:0007129            2.998713e-05                0.9999986          5
    ## 925   GO:0002430            4.767284e-05                0.9999997          3
    ##       numInCat                                           term ontology
    ## 1945       383      transmembrane signaling receptor activity       MF
    ## 8998       455                    signaling receptor activity       MF
    ## 12558      455                  molecular transducer activity       MF
    ## 464         76                      peptide receptor activity       MF
    ## 3106        17       homologous chromosome pairing at meiosis       BP
    ## 925          4 complement receptor mediated signaling pathway       BP

## Question 5

Using the Bottomly data, perform a differential expression analysis
using limma and treating strain as the outcome but adjusting for lane as
a factor. Then find genes significant at the 5% FDR rate using the
Benjamini Hochberg correction and perform the gene set analysis with
goseq following the protocol from the first 4 questions. How many of the
top 10 overrepresented categories are the same for the adjusted and
unadjusted analysis?

## Answer 5

First, we adjust the model for lane variable and get the list of
differentially expressed genes using the same

``` r
mod_adj = model.matrix(~pdata_bot$strain + as.factor(pdata_bot$lane.number))
fit_adj = lmFit(edata, mod_adj)
fit_adj = eBayes(fit_adj)
fit_BH_adj = topTable(fit_adj, number = dim(edata)[1], adjust.method = "BH", p.value = 0.05, sort.by="none")
```

    ## Removing intercept from test coefficients

Next, we create an integer list showing the DE (differentially expressed
at the 5% FDR level) status of each gene in the data set:

``` r
genes_adj = as.integer(rownames(edata) %in% rownames(fit_BH_adj))
names(genes_adj) = rownames(edata)
```

We then set up a weighting function for the genes:

``` r
pwf_adj = nullp(genes_adj,"mm9","ensGene")
```

    ## Loading mm9 length data...

![](Quiz-4_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Now, we perform the enrichment analysis:

``` r
GO.wall_adj = goseq(pwf_adj, "mm9", "ensGene")
```

    ## Fetching GO annotations...

    ## For 502 genes, we could not find any categories. These genes will be excluded.

    ## To force their use, please run with use_genes_without_cat=TRUE (see documentation).

    ## This was the default behavior for version 1.15.1 and earlier.

    ## Calculating the p-values...

    ## 'select()' returned 1:1 mapping between keys and columns

``` r
head(GO.wall_adj)
```

    ##         category over_represented_pvalue under_represented_pvalue numDEInCat
    ## 3106  GO:0007129            9.930757e-06                0.9999998          4
    ## 13571 GO:0070192            3.245292e-05                0.9999990          4
    ## 10253 GO:0045143            3.500942e-05                0.9999989          4
    ## 10247 GO:0045132            1.157886e-04                0.9999950          4
    ## 3104  GO:0007127            2.510831e-04                0.9999866          4
    ## 13435 GO:0061982            2.895429e-04                0.9999839          4
    ##       numInCat                                                   term ontology
    ## 3106        17               homologous chromosome pairing at meiosis       BP
    ## 13571       24 chromosome organization involved in meiotic cell cycle       BP
    ## 10253       25                      homologous chromosome segregation       BP
    ## 10247       35                         meiotic chromosome segregation       BP
    ## 3104        42                                              meiosis I       BP
    ## 13435       44                           meiosis I cell cycle process       BP

Finally we can find the number of common top 10 overrepresented
categories:

``` r
length(intersect(GO.wall$category[1:10], GO.wall_adj$category[1:10]))
```

    ## [1] 3

``` r
sessionInfo()
```

    ## R version 4.0.5 (2021-03-31)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 19042)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=English_United Kingdom.1252 
    ## [2] LC_CTYPE=English_United Kingdom.1252   
    ## [3] LC_MONETARY=English_United Kingdom.1252
    ## [4] LC_NUMERIC=C                           
    ## [5] LC_TIME=English_United Kingdom.1252    
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] org.Mm.eg.db_3.12.0         AnnotationDbi_1.52.0       
    ##  [3] DESeq2_1.30.1               SummarizedExperiment_1.20.0
    ##  [5] MatrixGenerics_1.2.1        matrixStats_0.59.0         
    ##  [7] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
    ##  [9] IRanges_2.24.1              S4Vectors_0.28.1           
    ## [11] limma_3.46.0                Biobase_2.50.0             
    ## [13] BiocGenerics_0.36.1         goseq_1.42.0               
    ## [15] geneLenDataBase_1.26.0      BiasedUrn_1.07             
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] nlme_3.1-152             bitops_1.0-7             bit64_4.0.5             
    ##  [4] RColorBrewer_1.1-2       progress_1.2.2           httr_1.4.2              
    ##  [7] tools_4.0.5              utf8_1.2.1               R6_2.5.0                
    ## [10] DBI_1.1.1                mgcv_1.8-36              colorspace_2.0-2        
    ## [13] tidyselect_1.1.1         prettyunits_1.1.1        bit_4.0.4               
    ## [16] curl_4.3.2               compiler_4.0.5           xml2_1.3.2              
    ## [19] DelayedArray_0.16.3      rtracklayer_1.49.5       scales_1.1.1            
    ## [22] genefilter_1.72.1        askpass_1.1              rappdirs_0.3.3          
    ## [25] stringr_1.4.0            digest_0.6.27            Rsamtools_2.6.0         
    ## [28] rmarkdown_2.9            XVector_0.30.0           pkgconfig_2.0.3         
    ## [31] htmltools_0.5.1.1        dbplyr_2.1.1             fastmap_1.1.0           
    ## [34] highr_0.9                rlang_0.4.11             RSQLite_2.2.7           
    ## [37] generics_0.1.0           BiocParallel_1.24.1      dplyr_1.0.7             
    ## [40] RCurl_1.98-1.3           magrittr_2.0.1           GO.db_3.12.1            
    ## [43] GenomeInfoDbData_1.2.4   Matrix_1.3-4             Rcpp_1.0.7              
    ## [46] munsell_0.5.0            fansi_0.5.0              lifecycle_1.0.0         
    ## [49] stringi_1.6.2            yaml_2.2.1               zlibbioc_1.36.0         
    ## [52] BiocFileCache_1.14.0     grid_4.0.5               blob_1.2.2              
    ## [55] crayon_1.4.1             lattice_0.20-44          Biostrings_2.58.0       
    ## [58] splines_4.0.5            GenomicFeatures_1.42.3   annotate_1.68.0         
    ## [61] hms_1.1.0                locfit_1.5-9.4           knitr_1.33              
    ## [64] pillar_1.6.1             geneplotter_1.68.0       biomaRt_2.46.3          
    ## [67] XML_3.99-0.6             glue_1.4.2               evaluate_0.14           
    ## [70] vctrs_0.3.8              gtable_0.3.0             openssl_1.4.4           
    ## [73] purrr_0.3.4              assertthat_0.2.1         cachem_1.0.5            
    ## [76] ggplot2_3.3.5            xfun_0.23                xtable_1.8-4            
    ## [79] survival_3.2-11          tibble_3.1.2             GenomicAlignments_1.26.0
    ## [82] memoise_2.0.0            ellipsis_0.3.2
