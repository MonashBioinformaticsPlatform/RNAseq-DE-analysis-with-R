# This file is generated from the corresponding .Rmd file



#### ---------------------------------
#### Learning Objectives {.objectives}
#### ---------------------------------

# =========================
# Install and load packages
# =========================

# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages. This
# requires internet connectivity and will take some time!
biocLite()


#### ---------------
#### NOTE {.callout}
#### ---------------

# ===================================
# Mapping reads to a reference genome
# ===================================

## -------------------
## Data pre-processing
## -------------------

# Delete all previously saved R objects
rm(list=ls())


# define shared directory for RNAseq data
RNAseqDATADIR <- "data/raw_data"
#list the fastq files in the raw data directory
dir(RNAseqDATADIR)


#### ---------------
#### NOTE {.callout}
#### ---------------

## -------
## Mapping
## -------

library(Rsubread)


#### ---------------
#### NOTE {.callout}
#### ---------------

#### 
####
#### 

source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")


# define the reference genome fasta file
REF_GENOME <- "hg19.fa"
# define the output directory for the Rsubread index
# (admin note: requires data/ref_data/download_hg19.sh to be run first)
RSUBREAD_INDEX_PATH <- "data/ref_data"
# define the basename for the index
RSUBREAD_INDEX_BASE <- "hg19"
# check what is in the reference directory
dir(RSUBREAD_INDEX_PATH)


#### ---------------
#### NOTE {.callout}
#### ---------------

#### 
####
#### 

# build the index
buildindex(basename=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), reference=REF_GENOME)


#### ---------------
#### NOTE {.callout}
#### ---------------

#### 
####
#### 

# list files in the raw data directory
dir(RNAseqDATADIR)
# define the fastq file with forward reads
inputfilefwd <- file.path(RNAseqDATADIR,"ERR420388_subsamp_1.fastq.gz")
# define the fastq file with reverse reads
inputfilervs <- file.path(RNAseqDATADIR,"ERR420388_subsamp_2.fastq.gz")

# run the align command to map the reads
align(index=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), readfile1=inputfilefwd, readfile2=inputfilervs, output_file="ERR420388.sam", output_format="SAM")


#### ----------------
#### NOTE: {.callout}
#### ----------------

#### 
####
#### 

# define the path to SAM file
outputsamfile <- "data/mapping/ERR420388.sam"
propmapped(outputsamfile)


# ============================
# Count reads for each feature
# ============================

# Getting read counts using the index previously built
mycounts<-featureCounts(outputsamfile, annot.ext=file.path(RSUBREAD_INDEX_PATH,"hg19.genes.gtf"), isGTFAnnotationFile=TRUE, isPairedEnd=TRUE)

# Checking your output object
summary(mycounts)
dim(mycounts$counts)
head(mycounts$annotation)
mycounts$targets
mycounts$stat



#### ----------------
#### NOTE: {.callout}
#### ----------------

#### 
####
#### 

MAPPINGDIR <- "data/mapping"
# load the counts previously calculated
load(file.path(MAPPINGDIR,"RawCounts.RData"))
# check the presence of read counts for the 8 libraries
summary(counts)


counts$targets


# print out counts table for every sample
write.table(counts$counts,file="~/raw_read_counts.txt",sep="\t", quote=F,append=F)


#### ----------------
#### NOTE: {.callout}
#### ----------------

#### 
####
#### 

# ============
# QC and stats
# ============

# define the experiment design file (tab separated text file is best)
EXPMT_DESIGN_FILE <- file.path(RNAseqDATADIR,'experiment_design.txt')
# read the experiment design file and save it into memory
experiment_design<-read.table(EXPMT_DESIGN_FILE,header=T,sep="\t")
#
# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- experiment_design$SampleID
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(counts$counts),]
# look at the design
experiment_design.ord
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$SampleID)
# create factors for future plotting
group<-factor(experiment_design.ord$tissue)
group
age<-factor(experiment_design.ord$age)
age


## --------------
## Basic QC plots
## --------------

#### ----------------
#### NOTE: {.callout}
#### ----------------

#### 
####
#### 

# density plot of raw read counts (log10)
png(file="~/Raw_read_counts_per_gene.density.png")
logcounts <- log(counts$counts[,1],10)
d <- density(logcounts)
plot(d,xlim=c(1,8),main="",ylim=c(0,.45),xlab="Raw read counts per gene (log10)", ylab="Density")
for (s in 2:length(samples)){
  logcounts <- log(counts$counts[,s],10)
  d <- density(logcounts)
  lines(d)
}
dev.off()


## box plots of raw read counts (log10)
png(file="~/Raw_read_counts_per_gene.boxplot.png")
logcounts <- log(counts$counts,10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:length(samples)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()


# select data for the 100 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts$counts[select,]

# heatmap with sample name on X-axis
png(file="~/High_expr_genes.heatmap.png")
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6))
dev.off()


# heatmap with condition group as labels
colnames(highexprgenes_counts)<- group
# plot
png(file="~/High_exprs_genes.heatmap.group.png")
heatmap(highexprgenes_counts, col = topo.colors(50), margin=c(10,6))
dev.off()


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
#
## ----------------------------
## Principal Component Analysis
## ----------------------------

# select data for the 1000 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:100]
highexprgenes_counts <- counts$counts[select,]
# annotate the data with condition group as labels
colnames(highexprgenes_counts)<- group
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
dim(data_for_PCA)


## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned


mds$eig


# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / sum(mds$eig)
# plot the PCA
png(file="~/PCA_PropExplainedVariance.png")
barplot(eig_pc,
     las=1,
     xlab="Dimensions",
     ylab="Proportion of explained variance (%)", y.axis=NULL,
     col="darkgrey")
dev.off()


## calculate MDS
mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis


#Samples representation
png(file="~/PCA_Dim1vsDim2.png")
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8)
dev.off()


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
# 
# 
#
# =======================
# Differential Expression
# =======================

# load required libraries
library(edgeR)


# get the expression counts from previous alignment step
mycounts <- counts$counts
dim(mycounts)
mycounts[1:5,1:3]

# filtering
#Keep genes with least 1 count-per-million reads (cpm) in at least 4 samples
isexpr <- rowSums(cpm(mycounts)>1) >= 4
table(isexpr)
mycounts <- mycounts[isexpr,]
genes <- rownames(mycounts)

dim(mycounts)


# load required libraries
library(limma)

# check your samples grouping
experiment_design.ord[colnames(mycounts),]$tissue == group

# create design matrix for limma
design <- model.matrix(~0+group)
# substitute "group" from the design column names
colnames(design)<- gsub("group","",colnames(design))
# check your design matrix
design

# calculate normalization factors between libraries
nf <- calcNormFactors(mycounts)
#
# normalise the read counts with 'voom' function
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
# extract the normalised read counts
counts.voom <- y$E

# save normalised expression data into output dir
write.table(counts.voom,file="~/counts.voom.txt",row.names=T,quote=F,sep="\t")

# fit linear model for each gene given a series of libraries
fit <- lmFit(y,design)
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(liver-brain,levels=design)
cont.matrix

# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)

# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
options(digits=3)

# check the output fit
dim(fit)


# set adjusted pvalue threshold and log fold change threshold
mypval=0.01
myfc=3

# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
mycoef="liver - brain"
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)

# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])

# get significant DE genes only (adjusted p-value < mypval)
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
dim(limma.res.pval)

# get significant DE genes with low adjusted p-value high fold change
limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
dim(limma.res.pval.FC)

# write limma output table for significant genes into a tab delimited file
filename = paste("~/DEgenes_",mycoef,"_pval",mypval,"_logFC",myfc,".txt",sep="")
write.table(limma.res.pval.FC,file=filename,row.names=T,quote=F,sep="\t")


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
# 
# 
# 
#
# ===============
# Gene Annotation
# ===============

# get the Ensembl annotation for human genome
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))

# get entrez gene IDs from limma output table
entrez_genes <- as.character(rownames(limma.res.pval.FC))
length(entrez_genes)

# interrogate the BioMart database to get gene symbol and description for these genes
detags.IDs <- getBM(
 filters= "entrezgene",
 attributes= c("entrezgene","hgnc_symbol","description"),
 values= entrez_genes,
 mart= mart
)

dim(detags.IDs)
head(detags.IDs)


# remove duplicates
detags.IDs.matrix<-detags.IDs[-which(duplicated(detags.IDs$entrezgene)),]
# select genes of interest only
rownames(detags.IDs.matrix)<-detags.IDs.matrix$entrezgene
entrez_genes.annot <- detags.IDs.matrix[as.character(entrez_genes),]
# join the two tables
rownames(limma.res.pval.FC) <- limma.res.pval.FC$ID
limma.res.pval.FC.annot <- cbind(entrez_genes.annot,limma.res.pval.FC)
# check the annotated table
head(limma.res.pval.FC.annot)


# ===================
# Gene Set Enrichment
# ===================

# load the library
library(GOstats)

# Define list of genes of interest (DE genes - EntrezGene IDs)
entrezgeneids <- as.character(rownames(limma.res.pval.FC))
length(entrezgeneids)
# Define the universe
universeids <- rownames(mycounts)
length(universeids)


# define the p-value cut off for the hypergeometric test
hgCutoff <- 0.05
params <- new("GOHyperGParams",annotation="org.Hs.eg",geneIds=entrezgeneids,universeGeneIds=universeids,ontology="BP",pvalueCutoff=hgCutoff,testDirection="over")
#  Run the test
hg <- hyperGTest(params)
# Check results
hg


## Get the p-values of the test
hg.pv <- pvalues(hg)
## Adjust p-values for multiple test (FDR)
hg.pv.fdr <- p.adjust(hg.pv,'fdr')
## select the GO terms with adjusted p-value less than the cut off
sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr < hgCutoff])
length(sigGO.ID)

# get table from HyperG test result
df <- summary(hg)
# keep only significant GO terms in the table
GOannot.table <- df[df[,1] %in% sigGO.ID,]
head(GOannot.table)


# Create text report of the significantly over-represented GO terms
write.table(GOannot.table,file="~/GOterms_OverRep_BP.txt",sep="\t",row.names=F)
# Create html report of all over-represented GO terms
htmlReport(hg, file="~/GOterms_OverRep_BP.html")


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
#### 
####
#### 

RDataFile <- "~/RNAseq_DE_analysis_with_R.RData"
save.image(RDataFile)


### ----------------------------------------------------
### Record package and version info with `sessionInfo()`
### ----------------------------------------------------

sessionInfo()

