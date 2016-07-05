---
output: html_document
---

# RNAseq DE analysis with R - Solutions

This document provides solutions to the exercises given in the 'RNAseq DE analysis with R' course.

## Clustering


```
## Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
## '/home/roxane/RNAseq_DE_analysis_with_R.RData', probable reason 'No such
## file or directory'
```

```
## Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

> ## Exercise: Heatmap  {.challenge}
> Produce a heatmap for the 50 most highly expressed genes and annotate the samples with with their age.
>
> * Subset the read counts object for the 50 most highly expressed genes
> * Annotate the samples in the subset with their age (check order with design!)
> * Plot a heatmap with this subset of data, scaling genes and ordering both genes and samples

Solution:


```r
# Subset the read counts object for the 30 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:50]
```

```
## Error in is.data.frame(x): object 'counts' not found
```

```r
highexprgenes_counts <- counts$counts[select,]
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

```r
# Annotate the samples in the subset with their age (check order with design!)
colnames(highexprgenes_counts)<- experiment_design.ord$age
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
head(highexprgenes_counts)
```

```
## Error in head(highexprgenes_counts): object 'highexprgenes_counts' not found
```

```r
# Plot a heatmap with this subset of data, scaling genes and ordering both genes and samples
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6))
```

```
## Error in heatmap(highexprgenes_counts, col = topo.colors(50), margin = c(10, : object 'highexprgenes_counts' not found
```


## Principal Component Analysis

> ## Exercise: PCA {.challenge}
> Produce a PCA plot from the read counts of the 500 most highly expressed genes and change the labels until you can identify the reason for the split between samples from the same tissue. 
>
> * Get the read counts for the 500 most highly expressed genes
> * Transpose this matrix of read counts
> * Check the number of dimensions explaining the variability in the dataset
> * Run the PCA with an appropriate number of components
> * Annotate the samples with their age \& re-run the PCA \& plot the main components
> * Annotate the samples with other clinical data \& re-run the PCA \& plot the main components until you can separate the samples within each tissue group


Solution:

```r
# select data for the 1000 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:500]
```

```
## Error in is.data.frame(x): object 'counts' not found
```

```r
highexprgenes_counts <- counts$counts[select,]
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

```r
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
```

```
## Error in t(highexprgenes_counts): object 'highexprgenes_counts' not found
```

```r
# Run the PCA with an appropriate number of components
mds <- cmdscale(dist(data_for_PCA))
```

```
## Error in as.matrix(x): object 'data_for_PCA' not found
```

```r
# Plot the PCA
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
```

```
## Error in plot(mds[, 1], -mds[, 2], type = "n", xlab = "Dimension 1", ylab = "Dimension 2", : object 'mds' not found
```

```r
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
```

```
## Error in text(mds[, 1], -mds[, 2], rownames(mds), cex = 0.8): object 'mds' not found
```

```r
# Annotate the samples with their age & re-run the PCA & plot the main components
rownames(mds) <- experiment_design.ord$age
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
```

```
## Error in plot(mds[, 1], -mds[, 2], type = "n", xlab = "Dimension 1", ylab = "Dimension 2", : object 'mds' not found
```

```r
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
```

```
## Error in text(mds[, 1], -mds[, 2], rownames(mds), cex = 0.8): object 'mds' not found
```

```r
# Annotate the samples with other clinical data & re-run the PCA & plot the main components until you can separate the samples within each tissue group
rownames(mds)<- experiment_design.ord$technical_replicate_group
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
```

```
## Error in plot(mds[, 1], -mds[, 2], type = "n", xlab = "Dimension 1", ylab = "Dimension 2", : object 'mds' not found
```

```r
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
```

```
## Error in text(mds[, 1], -mds[, 2], rownames(mds), cex = 0.8): object 'mds' not found
```


## Differential Expression

> ## Exercise: Limma {.challenge}
> Get the number of DE genes between technical group 1 and technical group 2 (all Brain samples) with adj pvalue<0.01.
>
> * Create a new design matrix for limma with the technical replicate groups
> * Re-normalise the read counts with 'voom' function with new design matrix
> * Fit a linear model on these normalised data
> * Make the contrast matrix corresponding to the new set of parameters
> * Fit the contrast matrix to the linear model
> * Compute moderated t-statistics of differential expression 
> * Get the output table for the 10 most significant DE genes for this comparison

Solution:


```r
library(limma)
```

```
## Loading required package: methods
```

```r
# Create a new design matrix for limma with the technical replicate groups
techgroup<-factor(experiment_design.ord$technical_replicate_group)
```

```
## Error in factor(experiment_design.ord$technical_replicate_group): object 'experiment_design.ord' not found
```

```r
design <- model.matrix(~0+techgroup)
```

```
## Error in eval(expr, envir, enclos): object 'techgroup' not found
```

```r
colnames(design)<- gsub("techgroup","",colnames(design))
```

```
## Error in is.data.frame(x): object 'design' not found
```

```r
design
```

```
## Error in eval(expr, envir, enclos): object 'design' not found
```

```r
# Re-normalise the read counts with 'voom' function with new design matrix
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
```

```
## Error in is(counts, "DGEList"): object 'mycounts' not found
```

```r
counts.voom <- y$E
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
# Fit a linear model on these normalised data
fit <- lmFit(y,design)
```

```
## Error in is(object, "list"): object 'y' not found
```

```r
# Make the contrast matrix corresponding to the new set of parameters
cont.matrix <- makeContrasts(group_2-group_1,levels=design)
```

```
## Error in is.factor(levels): object 'design' not found
```

```r
cont.matrix 
```

```
## Error in eval(expr, envir, enclos): object 'cont.matrix' not found
```

```r
# Fit the contrast matrix to the linear model
fit <- contrasts.fit(fit, cont.matrix)
```

```
## Error in NCOL(fit$coefficients): object 'fit' not found
```

```r
# Compute moderated t-statistics of differential expression 
fit <- eBayes(fit)
```

```
## Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, : object 'fit' not found
```

```r
options(digits=3)

# Get the output table for the 10 most significant DE genes for this comparison
dim(topTable(fit,coef="group_2 - group_1",p.val=0.01,n=Inf))
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
topTable(fit,coef="group_2 - group_1",p.val=0.01)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```


## Gene Ontology

> ## Exercise: GOstats {.challenge}
> Identify the GO terms in the Molecular Function domain that are over-represented (pvalue<0.01) in your list of DE genes.
>
> * Get your list of DE genes (Entrez Gene IDs)
> * Set the new parameters for the hypergeometric test
> * Run the test and adjust the pvalues in the output object
> * Identify the significant GO terms at pvalue 0.01

Solution:


```r
library(GOstats)
```

```
## Error in library(GOstats): there is no package called 'GOstats'
```

```r
# Get your list of DE genes (Entrez Gene IDs)
entrezgeneids <- as.character(rownames(limma.res.pval.FC))
```

```
## Error in rownames(limma.res.pval.FC): object 'limma.res.pval.FC' not found
```

```r
universeids <- rownames(mycounts)
```

```
## Error in rownames(mycounts): object 'mycounts' not found
```

```r
# Set the new parameters for the hypergeometric test
params <- new("GOHyperGParams",annotation="org.Hs.eg",geneIds=entrezgeneids,universeGeneIds=universeids,ontology="MF",pvalueCutoff=0.01,testDirection="over")
```

```
## Error in getClass(Class, where = topenv(parent.frame())): "GOHyperGParams" is not a defined class
```

```r
# Run the test and adjust the pvalues in the output object
hg <- hyperGTest(params)
```

```
## Error in eval(expr, envir, enclos): could not find function "hyperGTest"
```

```r
hg.pv <- pvalues(hg)
```

```
## Error in eval(expr, envir, enclos): could not find function "pvalues"
```

```r
hg.pv.fdr <- p.adjust(hg.pv,'fdr')
```

```
## Error in p.adjust(hg.pv, "fdr"): object 'hg.pv' not found
```

```r
# Identify the significant GO terms at pvalue 0.01
sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr < hgCutoff])
```

```
## Error in eval(expr, envir, enclos): object 'hg.pv.fdr' not found
```

```r
df <- summary(hg)
```

```
## Error in summary(hg): object 'hg' not found
```

```r
GOterms.sig <- df[df[,1] %in% sigGO.ID,"Term"]
```

```
## Error in df[, 1]: object of type 'closure' is not subsettable
```

```r
length(GOterms.sig )
```

```
## Error in eval(expr, envir, enclos): object 'GOterms.sig' not found
```

```r
head(GOterms.sig)
```

```
## Error in head(GOterms.sig): object 'GOterms.sig' not found
```


## R environment


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 14.04.1 LTS
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
## [1] methods   stats     graphics  grDevices utils     datasets  base     
## 
## other attached packages:
## [1] limma_3.16.7
## 
## loaded via a namespace (and not attached):
## [1] formatR_1.4   tools_3.2.0   knitr_1.13    stringr_0.6.2 evaluate_0.9
```
