# This file is generated from the corresponding .Rmd file



## ----------
## Clustering
## ----------

# Delete all previously saved R objects
rm(list=ls())

load("data/RNAseq_DE_analysis_with_R.RData")


##### -----------------
##### Exercise: Heatmap
##### -----------------
# Produce a heatmap for the 50 most highly expressed genes and annotate
# the samples with with their age.
# 
# * Subset the read counts object for the 50 most highly expressed genes
# * Annotate the samples in the subset with their age (check order with
# design!)
# * Plot a heatmap with this subset of data, scaling genes and ordering
# both genes and samples
# 
#
##### ---------
##### Solution:
##### ---------

# Subset the read counts object for the 30 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:50]
highexprgenes_counts <- counts$counts[select,]

# Annotate the samples in the subset with their age (check order with design!)
colnames(highexprgenes_counts)<- experiment_design.ord$age
head(highexprgenes_counts)

# Plot a heatmap with this subset of data, scaling genes and ordering both genes and samples
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6))


## ----------------------------
## Principal Component Analysis
## ----------------------------

##### -------------
##### Exercise: PCA
##### -------------
# Produce a PCA plot from the read counts of the 500 most highly
# expressed genes and change the labels until you can identify the
# reason for the split between samples from the same tissue.
# 
# * Get the read counts for the 500 most highly expressed genes
# * Transpose this matrix of read counts
# * Check the number of dimensions explaining the variability in the
# dataset
# * Run the PCA with an appropriate number of components
# * Annotate the samples with their age \& re-run the PCA \& plot the
# main components
# * Annotate the samples with other clinical data \& re-run the PCA \&
# plot the main components until you can separate the samples within
# each tissue group
# 
# 
#
##### ---------
##### Solution:
##### ---------

# select data for the 1000 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:500]
highexprgenes_counts <- counts$counts[select,]

# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)

# Run the PCA with an appropriate number of components
mds <- cmdscale(dist(data_for_PCA))

# Plot the PCA
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8)

# Annotate the samples with their age & re-run the PCA & plot the main components
rownames(mds) <- experiment_design.ord$age
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8)

# Annotate the samples with other clinical data & re-run the PCA & plot the main components until you can separate the samples within each tissue group
rownames(mds)<- experiment_design.ord$technical_replicate_group
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8)


## -----------------------
## Differential Expression
## -----------------------

##### ---------------
##### Exercise: Limma
##### ---------------
# Get the number of DE genes between technical group 1 and technical
# group 2 (all Brain samples) with adj pvalue<0.01.
# 
# * Create a new design matrix for limma with the technical replicate
# groups
# * Re-normalise the read counts with 'voom' function with new design
# matrix
# * Fit a linear model on these normalised data
# * Make the contrast matrix corresponding to the new set of parameters
# * Fit the contrast matrix to the linear model
# * Compute moderated t-statistics of differential expression
# * Get the output table for the 10 most significant DE genes for this
# comparison
# 
#
##### ---------
##### Solution:
##### ---------

library(limma)

# Create a new design matrix for limma with the technical replicate groups
techgroup<-factor(experiment_design.ord$technical_replicate_group)
design <- model.matrix(~0+techgroup)
colnames(design)<- gsub("techgroup","",colnames(design))
design

# Re-normalise the read counts with 'voom' function with new design matrix
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
counts.voom <- y$E

# Fit a linear model on these normalised data
fit <- lmFit(y,design)

# Make the contrast matrix corresponding to the new set of parameters
cont.matrix <- makeContrasts(group_2-group_1,levels=design)
cont.matrix

# Fit the contrast matrix to the linear model
fit <- contrasts.fit(fit, cont.matrix)

# Compute moderated t-statistics of differential expression
fit <- eBayes(fit)
options(digits=3)

# Get the output table for the 10 most significant DE genes for this comparison
dim(topTable(fit,coef="group_2 - group_1",p.val=0.01,n=Inf))
topTable(fit,coef="group_2 - group_1",p.val=0.01)


## -------------
## Gene Ontology
## -------------

##### -----------------
##### Exercise: GOstats
##### -----------------
# Identify the GO terms in the Molecular Function domain that are over-
# represented (pvalue<0.01) in your list of DE genes.
# 
# * Get your list of DE genes (Entrez Gene IDs)
# * Set the new parameters for the hypergeometric test
# * Run the test and adjust the pvalues in the output object
# * Identify the significant GO terms at pvalue 0.01
# 
# 
# 
#
#### ---------
#### Solution:
#### ---------

library(GOstats)
# Get your list of DE genes (Entrez Gene IDs)
entrezgeneids <- as.character(rownames(limma.res.pval.FC))
universeids <- rownames(mycounts)

# Set the new parameters for the hypergeometric test
params <- new("GOHyperGParams",annotation="org.Hs.eg",geneIds=entrezgeneids,universeGeneIds=universeids,ontology="MF",pvalueCutoff=0.01,testDirection="over")

# Run the test and adjust the pvalues in the output object
hg <- hyperGTest(params)
hg.pv <- pvalues(hg)
hg.pv.fdr <- p.adjust(hg.pv,'fdr')

# Identify the significant GO terms at pvalue 0.01
sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr < hgCutoff])
df <- summary(hg)
GOterms.sig <- df[df[,1] %in% sigGO.ID,"Term"]

length(GOterms.sig )
head(GOterms.sig)


#### -------------
#### R environment
#### -------------

sessionInfo()

