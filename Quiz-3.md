Coursera Statistical Genomics: Quiz 3
================

``` r
knitr::opts_chunk$set(cache = TRUE)
```

## Dependencies

This document depends on the following packages:

``` r
library(devtools)
library(Biobase)
library(SummarizedExperiment)
library(RSkittleBrewer)
library(broom)
library(limma)
library(snpStats)
library(plotrix)
library(genefilter)
library(DESeq2)
```

## Question 1

Load the example SNP data with the following code:

``` r
data(for.exercise)
use <- seq(1, ncol(snps.10), 10)
sub.10 <- snps.10[,use]
snpdata = sub.10@.Data
status = subject.support$cc
```

Fit a linear model and a logistic regression model to the data for the
3rd SNP. What are the coefficients for the SNP variable? How are they
interpreted? (Hint: Don’t forget to recode the 0 values to NA for the
SNP data)

## Answer 1

First we can visualise the data to gain an idea of what shape the data
form. Here I use the sizeplot function to represent the amount of
overlapping points on the scatterpoint graph:

``` r
sizeplot(snpdata[,3], status, col =2, scale =.4, pch = 19, ylim = c(-0.1, 1.1), ylab = "Status", xlab = "SNP")
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Next, we can fit a linear model to the data:

``` r
lm3 = lm(snpdata[,3] ~ status)
tidy(lm3)
```

    ## # A tibble: 2 x 5
    ##   term        estimate std.error statistic p.value
    ##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 (Intercept)   1.12      0.0159    70.3     0    
    ## 2 status       -0.0180    0.0225    -0.799   0.424

Before we can fit a logistic regression model ‘0’ values must be changed
to ‘NA’ values:

``` r
snp3 = as.numeric(snpdata[,3])
snp3[snp3==0] = NA
```

Now we can fit the model:

``` r
glm3 = glm(status ~ snp3,family="binomial")
tidy(glm3)
```

    ## # A tibble: 2 x 5
    ##   term        estimate std.error statistic p.value
    ##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 (Intercept)    0.177     0.220     0.806   0.420
    ## 2 snp3          -0.158     0.188    -0.841   0.400

Both models are fit on the additive scale. So in the linear model case,
the coefficient is the decrease in probability associated with each
additional copy of the minor allele. In the logistic regression case, it
is the decrease in the log odds ratio associated with each additional
copy of the minor allele.

## Question 2

In the previous question why might the choice of logistic regression be
better than the choice of linear regression?

## Answer 2

If you included more variables it would be possible to get negative
estimates for the probability of being a case from the linear model, but
this would be prevented with the logistic regression model.

## Question 3

Fit a logistic regression model on a recessive (need 2 copies of minor
allele to confer risk) and additive scale for the 10th SNP. Make a table
of the fitted values versus the case/control status. Does one model fit
better than the other?

## Answer 3

First, we fit the recessive model.

``` r
snp10 = as.numeric(snpdata[,10])
snp10_rec = (snp10 == 3)
glm10_rec = glm(status ~ snp10_rec,family="binomial")
tidy(glm10_rec)
```

    ## # A tibble: 2 x 5
    ##   term          estimate std.error statistic p.value
    ##   <chr>            <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 (Intercept)    0.00922    0.0679     0.136   0.892
    ## 2 snp10_recTRUE -0.0698     0.187     -0.374   0.709

Then, the additive model.

``` r
snp10[snp10==0] = NA
glm10 = glm(status ~ snp10,family="binomial")
tidy(glm10)
```

    ## # A tibble: 2 x 5
    ##   term        estimate std.error statistic p.value
    ##   <chr>          <dbl>     <dbl>     <dbl>   <dbl>
    ## 1 (Intercept) -0.00751    0.174    -0.0433   0.965
    ## 2 snp10        0.00201    0.0933    0.0215   0.983

Next we can plot the fitted values against their case-control status (0
= Control, 1 = Case)

``` r
table(glm10_rec$fitted.values, status)
```

    ##                    status
    ##                       0   1
    ##   0.484848484848519  68  64
    ##   0.502304147465438 432 436

As NA/0 values have been excluded automatically, we must also exclude
them from the status object before comparing them against fitted values.

``` r
status_noNA = status[-which(snp10 %in% NA)]
table(glm10$fitted.values, status_noNA)
```

    ##                    status_noNA
    ##                       0   1
    ##   0.498624476867283 202 197
    ##   0.499127261301324 227 234
    ##   0.499630047500349  68  64

``` r
sizeplot(glm10_rec$fitted.values, status, col =2, scale =.5, pch = 19, xlim = c(0.48, 0.51), ylim = c(-0.2, 1.2), main = "Recessive Model", ylab = "Status", xlab = "Fitted Values: recessive model")
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-11-1.png)<!-- -->

``` r
sizeplot(glm10$fitted.values, status_noNA, col =2, scale =.5, pch = 19, xlim = c(0.4985, 0.5), ylim = c(-0.2, 1.2), main = "Additive Model", ylab = "Status", xlab = "Fitted Values")
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-11-2.png)<!-- -->

The recessive model fits much better since it appears that once you
aggregate the heterozygotes and homozygous minor alleles, there is a
bigger difference in the proportion of cases and controls.

## Question 4

Fit an additive logistic regression model to each SNP. What is the
average effect size? What is the max? What is the minimum?

## Answer 4

We can use a for loop to iterate through each gene, fitting the data and
extracting the effect size (z value) and storing it in a matrix.

``` r
z_score=matrix(NA,2851,1)

for (i in 1:2851){
    
    snp=as.numeric(snpdata[,i])
    
    snp[snp==0]=NA
    
    z_score[i,]=coef(summary(glm(status~snp,family="binomial")))[2,"z value"]
    
}
```

Next, compute the mean and the minimum and maximum value in the matrix.

``` r
mean(z_score)
```

    ## [1] 0.007155377

``` r
min(z_score)
```

    ## [1] -4.251469

``` r
max(z_score)
```

    ## [1] 3.900891

## Question 5

Fit an additive logistic regression model to each SNP and square the
coefficients. What is the correlation with the results from using
snp.rhs.tests and chi.squared? Why does this make sense?

## Answer 5

First we square the z\_score values we obtained in Question 4

``` r
z_score_squ = z_score^2
```

Then we perform snp.rhs.tests on the snp data and obtain the chi squared
values.Then find the correlation between the two.

``` r
glm_all = snp.rhs.tests(status ~ 1,snp.data=sub.10)

chi_values <- chi.squared(glm_all)

cor(z_score_squ, chi_values)
```

    ##           [,1]
    ## [1,] 0.9992946

They are both testing for the same association using the same additive
regression model on the logistic scale but using slightly different
tests, hence the correlation is not exactly 1.

## Question 6

Load the Montgomery and Pickrell eSet:

``` r
con =url("http://bowtie-bio.sourceforge.net/recount/ExpressionSets/montpick_eset.RData")
load(file=con)
close(con)
mp = montpick.eset
pdata=pData(mp)
edata=as.data.frame(exprs(mp))
fdata = fData(mp)
```

Do the log2(data + 1) transform and fit calculate F-statistics for the
difference between studies/populations using genefilter:rowFtests and
using genefilter:rowttests. Do you get the same statistic? Do you get
the same p-value?

## Answer 6

First we trandform the data:

``` r
edata_log2 = log2(edata + 1)
```

Next, we can perform t-tests and F-tests for each row of the log2
transformed edata with respect to the population variable:

``` r
tstats_obj = rowttests(as.matrix(edata_log2),as.factor(pdata$population))
fstats_obj = rowFtests(as.matrix(edata_log2),as.factor(pdata$population))
```

Next we plot the t statistics vs the F statistics. If they are identical
tests we should expect to see a straight line with a positive gardient
of 1.

``` r
plot(tstats_obj$statistic, fstats_obj$statistic)
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-19-1.png)<!-- -->

We do not see a straight line but a polynomial confirming that, in this
case, the F statistic is a transform of the t statistic.

Next, we can plot the P-values of each test:

``` r
plot(tstats_obj$p.value, fstats_obj$p.value)
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-20-1.png)<!-- -->

We get the same p-value but different statistics between tests. This is
because the F-statistic and t-statistic test the same thing when doing a
two group test and one is a transform of the other.

## Question 7

Using the following rowMeans &gt; 100 Montgomery and Pickrell eSet:

``` r
edata = edata[rowMeans(edata) > 100,]
```

First test for differences between the studies using the DESeq2 package
using the DESeq function. Then do the log2(data + 1) transform and do
the test for differences between studies using the limma package and the
lmFit, ebayes and topTable functions. What is the correlation in the
statistics between the two analyses? Are there more differences for the
large statistics or the small statistics (hint: Make an MA-plot).

## Answer 7

First, a regression is performed using the DESeq2 package:

``` r
de = DESeqDataSetFromMatrix(countData = edata, colData = pdata, design =  ~population)
glm_all_de = DESeq(de)
```

    ## estimating size factors

    ## estimating dispersions

    ## gene-wise dispersion estimates

    ## mean-dispersion relationship

    ## final dispersion estimates

    ## fitting model and testing

    ## -- replacing outliers and refitting for 3 genes
    ## -- DESeq argument 'minReplicatesForReplace' = 7 
    ## -- original counts are preserved in counts(dds)

    ## estimating dispersions

    ## fitting model and testing

``` r
result_de = results(glm_all_de)
```

Then, after performing a log2 transform, we perform linear regressions:

``` r
edata_log2 = log2(edata + 1)
mod_study = model.matrix(~ as.factor(pdata$population))
fit_limma_study = lmFit(edata_log2, mod_study)
ebayes_limma_study = eBayes(fit_limma_study)
result_lim = topTable(ebayes_limma_study, number=nrow(edata_log2))
```

    ## Removing intercept from test coefficients

The results from each regression are confirmed to be in a matching order
so they can be compared and their correlation is calculated.

``` r
result_lim=result_lim[match(rownames(result_de),rownames(result_lim)),]
cor(result_de$stat, result_lim$t)
```

    ## [1] 0.9278568

Next, we can make some MA-plots to see the distribution of logFC against
the statistics.

``` r
par(mfcol=c(1,2))
plot(log2(result_de$baseMean+1), result_de$log2FoldChange, ylim=c(-5,5), pch=16)
plot(result_lim$AveExpr, result_lim$logFC, ylim=c(-5,5), pch=16)
```

![](Quiz-3_files/figure-gfm/unnamed-chunk-25-1.png)<!-- -->

As shown in the two MA-plots above, there are more differences in logFC
for smaller statistics.

## Question 8

Apply the Benjamni-Hochberg correction to the P-values from the two
previous analyses. How many results are statistically significant at an
FDR of 0.05 in each analysis?

## Answer 8

Benjamni-Hochberg corrections are applied to the P-values of each
experiment and then the number of the adjusted P-values below 0.05 are
counted:

``` r
deseq2_fdr = p.adjust(result_de$pvalue,"BH")
deseq2_fdr = deseq2_fdr[!is.na(deseq2_fdr)]

limma_fdr = p.adjust(result_lim$P.Value,"BH")
limma_fdr = limma_fdr[!is.na(limma_fdr)]

print(length(deseq2_fdr[deseq2_fdr < 0.05]))
```

    ## [1] 1995

``` r
print(length(limma_fdr[limma_fdr < 0.05]))
```

    ## [1] 2807

## Question 9

Is the number of significant differences surprising for the analysis
comparing studies from Question 8? Why or why not?

## Answer 9

One one hand, it is surprising as there is a large fraction of the genes
that are significantly different, but it isn’t that surprising when one
considers that comparisons are being made between measurements from very
different batches.

## Question 10

Suppose you observed the following P-values from the comparison of
differences between studies. Why might you be suspicious of the
analysis?

<img src="images/pvals.png" width="90%" />

## Answer 10

The p-values should have a spike near zero (the significant results) and
be flat to the right hand side (the null results) so the distribution
pushed toward one suggests something went wrong.

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
    ##  [1] DESeq2_1.30.1               genefilter_1.72.1          
    ##  [3] plotrix_3.8-1               snpStats_1.40.0            
    ##  [5] Matrix_1.3-4                survival_3.2-11            
    ##  [7] limma_3.46.0                broom_0.7.6                
    ##  [9] RSkittleBrewer_1.1          SummarizedExperiment_1.20.0
    ## [11] GenomicRanges_1.42.0        GenomeInfoDb_1.26.7        
    ## [13] IRanges_2.24.1              S4Vectors_0.28.1           
    ## [15] MatrixGenerics_1.2.1        matrixStats_0.59.0         
    ## [17] Biobase_2.50.0              BiocGenerics_0.36.1        
    ## [19] devtools_2.4.2              usethis_2.0.1              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] bitops_1.0-7           fs_1.5.0               bit64_4.0.5           
    ##  [4] RColorBrewer_1.1-2     httr_1.4.2             rprojroot_2.0.2       
    ##  [7] tools_4.0.5            backports_1.2.1        utf8_1.2.1            
    ## [10] R6_2.5.0               DBI_1.1.1              colorspace_2.0-2      
    ## [13] withr_2.4.2            tidyselect_1.1.1       prettyunits_1.1.1     
    ## [16] processx_3.5.2         bit_4.0.4              compiler_4.0.5        
    ## [19] cli_3.0.0              desc_1.3.0             DelayedArray_0.16.3   
    ## [22] scales_1.1.1           callr_3.7.0            stringr_1.4.0         
    ## [25] digest_0.6.27          rmarkdown_2.9          XVector_0.30.0        
    ## [28] pkgconfig_2.0.3        htmltools_0.5.1.1      sessioninfo_1.1.1     
    ## [31] highr_0.9              fastmap_1.1.0          rlang_0.4.11          
    ## [34] rstudioapi_0.13        RSQLite_2.2.7          generics_0.1.0        
    ## [37] BiocParallel_1.24.1    dplyr_1.0.6            RCurl_1.98-1.3        
    ## [40] magrittr_2.0.1         GenomeInfoDbData_1.2.4 Rcpp_1.0.6            
    ## [43] munsell_0.5.0          fansi_0.5.0            lifecycle_1.0.0       
    ## [46] stringi_1.6.2          yaml_2.2.1             zlibbioc_1.36.0       
    ## [49] pkgbuild_1.2.0         grid_4.0.5             blob_1.2.1            
    ## [52] crayon_1.4.1           lattice_0.20-44        splines_4.0.5         
    ## [55] annotate_1.68.0        locfit_1.5-9.4         knitr_1.33            
    ## [58] ps_1.6.0               pillar_1.6.1           geneplotter_1.68.0    
    ## [61] codetools_0.2-18       pkgload_1.2.1          XML_3.99-0.6          
    ## [64] glue_1.4.2             evaluate_0.14          remotes_2.4.0         
    ## [67] vctrs_0.3.8            testthat_3.0.2         gtable_0.3.0          
    ## [70] purrr_0.3.4            tidyr_1.1.3            assertthat_0.2.1      
    ## [73] cachem_1.0.5           ggplot2_3.3.5          xfun_0.23             
    ## [76] xtable_1.8-4           tibble_3.1.2           AnnotationDbi_1.52.0  
    ## [79] memoise_2.0.0          ellipsis_0.3.2
