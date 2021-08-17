Coursera Statistical Genomics: Quiz 2
================

``` r
knitr::opts_chunk$set(cache = TRUE)
```

## Dependencies

This document depends on the following packages:

``` r
library(devtools)
```

    ## Loading required package: usethis

``` r
library(Biobase)
```

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

``` r
library(SummarizedExperiment)
```

    ## Loading required package: MatrixGenerics

    ## Loading required package: matrixStats

    ## 
    ## Attaching package: 'matrixStats'

    ## The following objects are masked from 'package:Biobase':
    ## 
    ##     anyMissing, rowMedians

    ## 
    ## Attaching package: 'MatrixGenerics'

    ## The following objects are masked from 'package:matrixStats':
    ## 
    ##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
    ##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
    ##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
    ##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
    ##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
    ##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
    ##     colWeightedMeans, colWeightedMedians, colWeightedSds,
    ##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
    ##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
    ##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
    ##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
    ##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
    ##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
    ##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
    ##     rowWeightedSds, rowWeightedVars

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     rowMedians

    ## Loading required package: GenomicRanges

    ## Loading required package: stats4

    ## Loading required package: S4Vectors

    ## 
    ## Attaching package: 'S4Vectors'

    ## The following object is masked from 'package:base':
    ## 
    ##     expand.grid

    ## Loading required package: IRanges

    ## 
    ## Attaching package: 'IRanges'

    ## The following object is masked from 'package:grDevices':
    ## 
    ##     windows

    ## Loading required package: GenomeInfoDb

``` r
library(RSkittleBrewer)
library(broom)
library(limma)
```

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA

## Question 1

Load the Montgomery and Pickrell eSet (code supplied by Coursera):

``` r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
```

What percentage of variation is explained by the 1st principal component
in the data set if you:

1.  Do no transformations?
2.  log2(data + 1) transform?
3.  log2(data + 1) transform and subtract row means?

# Answer 1

First, the appropriate data transforms are applied to the dataset:

``` r
edata_log2 = log2(edata + 1)
edata_RM = edata_log2 - rowMeans(edata_log2)
```

Then, the Singular-Value Decompositions (SVD) are calculated:

``` r
svd1 = svd(edata)
svd2 = svd(edata_log2)
svd3 = svd(edata_RM)
```

The percent variance explained are shown:

``` r
plot(svd1$d,ylab="Singular value",col=1)
```

![](Week-2_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
plot(svd2$d,ylab="Singular value",col=2)
```

![](Week-2_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
plot(svd3$d,ylab="Singular value",col=3)
```

![](Week-2_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

From the plots shown above, we can expect PC1 of svd2 to acount for the
highest variance, followed by svd1 and, finally, svd3.

Next, we work out the exact variance explained by each PC1:

``` r
(svd1$d^2/sum(svd1$d^2))[1]
```

    ## [1] 0.8873421

``` r
(svd2$d^2/sum(svd2$d^2))[1]
```

    ## [1] 0.9737781

``` r
(svd3$d^2/sum(svd3$d^2))[1]
```

    ## [1] 0.3463729

## Question 2

Using the Montgomery and Pickerell eSet, perform the log2(data + 1)
transform and subtract row means from the samples. Set the seed to 333
and use k-means to cluster the samples into two clusters. Use svd to
calculate the singular vectors. What is the correlation between the
first singular vector and the sample clustering indicator?

## Answer 2

as the log2 and row means transform has already been applied
(edata\_RM), first step is to set the seed and run the kmeans with 2
centres:

``` r
set.seed(333)
kmeans1 = kmeans(t(edata_RM), centers = 2)
```

Next we find the correlation between the first singular vector and the
sample clustering indicator:

``` r
cor(svd3$v[,1], kmeans1$cluster)
```

    ## [1] -0.8678247

## Question 3

Question 3

Load the Bodymap data with the following command (code below supplied by
Coursera):

``` r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/bodymap_eset.RData")
load(file=con)
close(con)
bm = bodymap.eset
edata_bm = exprs(bm)
pdata_bm=pData(bm)
```

Fit a linear model relating the first gene’s counts to the number of
technical replicates, treating the number of replicates as a factor.
Plot the data for this gene versus the covariate. Can you think of why
this model might not fit well?

## Answer 3

First we fit a linear model to the first gene’s counts and the number of
technical replicates per sample.

``` r
edata_bm = as.matrix(edata_bm)
lm1 = lm(edata_bm[1,] ~ as.factor(pdata_bm$num.tech.reps))
```

Then plot the data points on a gene 1 counts vs technical reps graph and
overlay the linear regression line to visualize how the model fits to
the data.

``` r
plot(pdata_bm$num.tech.reps, edata_bm[1,], col = 1)
abline(lm1$coefficients[1], lm1$coefficients[2], col = 2, lwd = 3)
```

![](Week-2_files/figure-gfm/unnamed-chunk-12-1.png)<!-- --> (There may
be different numbers of counts for different numbers of technical
replicates.) The plot shows that the data are highly skewed to the
right.

## Question 4

Using the Bodymap data, fit a linear model relating the first gene’s
counts to the age of the person and the sex of the samples. What is the
value and interpretation of the coefficient for age?

## Answer 4

Simply, run a linear regression model for the first gene in edata\_bm
relating to age and sex of the samples using the lm() function:

``` r
lm2 = lm(edata_bm[1,] ~ pdata_bm$age + pdata_bm$gender)
lm2
```

    ## 
    ## Call:
    ## lm(formula = edata_bm[1, ] ~ pdata_bm$age + pdata_bm$gender)
    ## 
    ## Coefficients:
    ##      (Intercept)      pdata_bm$age  pdata_bm$genderM  
    ##          2331.58            -23.91           -207.26

The coefficient for age is -23.91, meaning that for each year of age the
counts go down by 23.91 assuming gender is fixed.

## Question 5

Using the Montgomery and Pickrell eSet, perform the log2(data + 1)
transform. Then fit a regression model to each sample using population
as the outcome. Do this using the lm.fit function. What is the dimension
of the residual matrix, the effects matrix and the coefficients matrix?

## Answer 5

First create a model matrix using population, then fit a linear model to
the matrix and the log2 + 1 transformed expression data:

``` r
mod = model.matrix(~ pdata$population)
fit = lm.fit(mod, t(edata))
```

Finally, show the dimensions of the residual matrix, the effects matrix
and the coefficients matrix:

``` r
dim(fit$residuals)
```

    ## [1]   129 52580

``` r
dim(fit$effects)
```

    ## [1]   129 52580

``` r
dim(fit$coefficients)
```

    ## [1]     2 52580

## Question 6

Using the Montgomery and Pickrell eSet, perform the log2(data + 1)
transform. Then fit a regression model to each sample using population
as the outcome. Do this using the lm.fit function. What is the effects
matrix?

## Answer 6

Using the fit object created in question 5:

``` r
fit$effects[1,1:5]
```

    ## ENSG00000000003 ENSG00000000005 ENSG00000000419 ENSG00000000457 ENSG00000000460 
    ##      -1.2326313      -0.1760902    -545.9676070    -642.2008911    -123.7913974

By reading the documentation for lm.fit, we can infer that the effects
matrix shows the estimated fitted values for all samples for each gene,
with the values for each gene stored in the columns of the matrix.

## Question 7

Using the Bodymap dataset, fit many regression models to the expression
data where age is the outcome variable using the lmFit function from the
limma package (hint: you may have to subset the expression data to the
samples without missing values of age to get the model to fit). What is
the coefficient for age for the 1,000th gene? Make a plot of the data
and fitted values for this gene. Does the model fit well?

## Answer 7

First, check for NA values and remove the samples containing these
missing data from edata and pdata by their indices:

``` r
table(pdata_bm$age,useNA="ifany")
```

    ## 
    ##   19   29   37   47   58   60   65   68   73   77   86 <NA> 
    ##    1    1    1    1    1    3    1    1    2    3    1    3

``` r
excIndex = which(is.na(pdata_bm$age))
edata_bm_cl = edata_bm[,-excIndex]
pdata_bm_cl = pdata_bm[-excIndex,]
```

Now we check the dimensions of the new edata and pdata to make sure that
they match:

``` r
dim(edata_bm_cl)
```

    ## [1] 52580    16

``` r
dim(pdata_bm_cl)
```

    ## [1] 16  6

First create a model matrix using age as a variable:

``` r
mod2 = model.matrix(~ pdata_bm_cl$age)
```

Fit a regression to each gene using age as the outcome variable using
the limma package function lmFit() and identify the coefficient for age
of the 1000th gene:

``` r
fit_limma = lmFit(edata_bm_cl, mod2)
fit_limma$coefficients[1000,]
```

    ##     (Intercept) pdata_bm_cl$age 
    ##      2469.87375       -27.61178

Now we can plot the data as well as the fitted line:

``` r
plot(pdata_bm_cl$age, edata_bm_cl[1000,], col =2)
abline(fit_limma[1000,], col = 1, lwd = 3)
```

![](Week-2_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

The model clearly does not fit well due to two large outlying values
while the rest of the values lie near zero.

## Question 8

Using the Bodymap data, fit many regression models to the expression
data where age is the outcome variable and tissue.type is an adjustment
variable using the lmFit function from the limma package (hint: you may
have to subset the expression data to the samples without missing values
of age to get the model to fit). What is wrong with this model?

## Answer 8

``` r
mod3 = model.matrix(~ pdata_bm_cl$age + as.numeric(pdata_bm_cl$tissue.type))
fit_limma2 = lmFit(edata_bm_cl, mod3)
```

## Session Info

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
    ##  [1] limma_3.46.0                broom_0.7.6                
    ##  [3] RSkittleBrewer_1.1          SummarizedExperiment_1.20.0
    ##  [5] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
    ##  [7] IRanges_2.24.1              S4Vectors_0.28.1           
    ##  [9] MatrixGenerics_1.2.1        matrixStats_0.59.0         
    ## [11] Biobase_2.50.0              BiocGenerics_0.36.1        
    ## [13] devtools_2.4.2              usethis_2.0.1              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] lattice_0.20-44        tidyr_1.1.3            prettyunits_1.1.1     
    ##  [4] ps_1.6.0               assertthat_0.2.1       rprojroot_2.0.2       
    ##  [7] digest_0.6.27          utf8_1.2.1             R6_2.5.0              
    ## [10] backports_1.2.1        evaluate_0.14          highr_0.9             
    ## [13] pillar_1.6.1           zlibbioc_1.36.0        rlang_0.4.11          
    ## [16] callr_3.7.0            Matrix_1.3-4           rmarkdown_2.9         
    ## [19] desc_1.3.0             stringr_1.4.0          RCurl_1.98-1.3        
    ## [22] DelayedArray_0.16.3    compiler_4.0.5         xfun_0.23             
    ## [25] pkgconfig_2.0.3        pkgbuild_1.2.0         htmltools_0.5.1.1     
    ## [28] tidyselect_1.1.1       tibble_3.1.2           GenomeInfoDbData_1.2.4
    ## [31] codetools_0.2-18       fansi_0.5.0            crayon_1.4.1          
    ## [34] dplyr_1.0.6            withr_2.4.2            bitops_1.0-7          
    ## [37] grid_4.0.5             lifecycle_1.0.0        DBI_1.1.1             
    ## [40] magrittr_2.0.1         cli_3.0.0              stringi_1.6.2         
    ## [43] cachem_1.0.5           XVector_0.30.0         fs_1.5.0              
    ## [46] remotes_2.4.0          testthat_3.0.2         ellipsis_0.3.2        
    ## [49] vctrs_0.3.8            generics_0.1.0         tools_4.0.5           
    ## [52] glue_1.4.2             purrr_0.3.4            processx_3.5.2        
    ## [55] pkgload_1.2.1          fastmap_1.1.0          yaml_2.2.1            
    ## [58] sessioninfo_1.1.1      memoise_2.0.0          knitr_1.33
