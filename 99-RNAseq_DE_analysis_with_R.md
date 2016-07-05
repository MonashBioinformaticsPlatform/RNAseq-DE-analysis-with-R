---
output: html_document
---

# RNA-seq: differential gene expression analysis  

> ## Learning Objectives {.objectives}
> This course is an introduction to differential expression analysis from RNAseq data. It will take you from the raw fastq files all the way to the list of differentially expressed genes, via the mapping of the reads to a reference genome and statistical analysis using the limma package.
>
> The Monash Bioinformatics Platform thanks and acknowledges QFAB (<http://www.qfab.org>) for original materials in this section.



### Install and load packages

Most generic R packages are hosted on the Comprehensive R Archive Network (CRAN, <http://cran.us.r-project.org/>). To install one of these packages, you would use `install.packages("packagename")`. You only need to install a package once, then load it each time using `library(packagename)`. 

Bioconductor packages work a bit differently, and are not hosted on CRAN. Go to <http://bioconductor.org/> to learn more about the Bioconductor project. To use any Bioconductor package, you'll need a few "core" Bioconductor packages. Run the following commands to (1) download the installer script, and (2) install some core Bioconductor packages. You'll need internet connectivity to do this, and it'll take a few minutes, but it only needs to be done once.


```r
# Download the installer script
source("http://bioconductor.org/biocLite.R")

# biocLite() is the bioconductor installer function. Run it without any
# arguments to install the core packages or update any installed packages. This
# requires internet connectivity and will take some time!
biocLite()
```


> ## NOTE {.callout}
>
> All packages required for this course have already been installed on our training server so you won't need to run this part today.


### Data pre-processing  

To start with, delete all previously saved R objects and define your working directory for the RNAseq data analysis.


```r
# Delete all previously saved R objects
rm(list=ls())
```

The data considered for the RNAseq part of the course have been downloaded from ArrayExpress (<http://www.ebi.ac.uk/arrayexpress>) and correspond to 8 RNA sequencing libraries from Human brain and liver.

Raw sequencing data are usually available in FASTQ format which is a well defined text-based format for storing both biological sequences (usually nucleotide sequences) and their corresponding quality scores. The raw data from this study have been downloaded (8Gb / fastq file) into the shared directory "~/data/RNAseq/raw_data".



```r
# define shared directory for RNAseq data
RNAseqDATADIR <- "/mnt/RNAseqCourse/raw_data"
#list the fastq files in the raw data directory
dir(RNAseqDATADIR)
```

```
## character(0)
```


The first step in a RNAseq analysis is to run a quick quality check on your data, this will give you an idea of the quality of your raw data in terms of number of reads per library, read length, average quality score along the reads, GC content, sequence duplication level, adaptors that might have not been removed correctly from the data etc.

The fastQC tool is quick and easy to run and can be downloaded from here: <http://www.bioinformatics.babraham.ac.uk/projects/fastqc/>.

To ensure highest quality of the sequences for subsequent mapping and differential expression analysis steps, the reads can also be trimmed using the Trimmomatic tool (Lohse et al. 2012, <http://www.usadellab.org/cms/?page=trimmomatic>).

 
> ## NOTE {.callout}
>
> For the scope of this course we will focus on the R-based steps and will assume that the data are fit for purpose.
>
> Also the next two sections (Mapping and Read counts) will be performed on small subsets of those files due to time and I/O limitations. 



### Mapping reads to a reference genome

Once the reads have been quality checked and trimmed, the next step is to map the reads to the reference genome (in our case the human genome "hg19"). This can be done with the Bioconductor package *Rsubread* (Y et al. 2013).


```r
library(Rsubread)
```

```
## Error in library(Rsubread): there is no package called 'Rsubread'
```


> ## NOTE {.callout}
> 
> Please remember that you can also perform those steps using command line tools (such as BWA and featureCounts) as described in the Unix shell course.


If the package is not already installed on your system, you can install it by typing:

```r
source("http://bioconductor.org/biocLite.R")
biocLite("Rsubread")
```


Rsubread provides reference genome indices for the most common organisms: human and mouse. If you are working with a different organism you can build your own index using the *buildindex* command.


```r
# define the reference genome fasta file
REF_GENOME <- "hg19.fa"
# define the output directory for the Rsubread index
RSUBREAD_INDEX_PATH <- "/mnt/RNAseqCourse/ref_data"
# define the basename for the index
RSUBREAD_INDEX_BASE <- "hg19"
# check what is in the reference directory
dir(RSUBREAD_INDEX_PATH)
```


> ## NOTE {.callout}
>
> For the purpose of this course the index has already been built and the mapping will be done on small subset of reads since this step can take up to a couple of hours per library. 
> **Please do not run the command below**.



```r
# build the index
buildindex(basename=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), reference=REF_GENOME)
```



Once the Rsubread index has been created you can map your reads to the genome by running the *align* command.

The code below will be used to map the reads for a specific library against the genome for which the index has been built.



```r
# list files in the raw data directory
dir(RNAseqDATADIR)
# define the fastq file with forward reads
inputfilefwd <- file.path(RNAseqDATADIR,"ERR420388_subsamp_1.fastq.gz")
# define the fastq file with reverse reads
inputfilervs <- file.path(RNAseqDATADIR,"ERR420388_subsamp_2.fastq.gz")

# run the align command to map the reads
align(index=file.path(RSUBREAD_INDEX_PATH,RSUBREAD_INDEX_BASE), readfile1=inputfilefwd, readfile2=inputfilervs, output_file="ERR420388.sam", output_format="SAM")
```


> ## NOTE: {.callout}
>
> The nthreads parameter can be used in the *align* command to speed up the process and run the alignment using several CPUs in parallel.


The function *propmapped* returns the proportion of mapped reads in the output SAM file: total number of input reads, number of mapped reads and proportion of mapped reads. Let's have a look at a small SAM file as an example:


```r
# define the path to SAM file  
outputsamfile <- "/mnt/RNAseqCourse/mapping/ERR420388.sam"
propmapped(outputsamfile)
```


### Count reads for each feature 

Rsubread provides a read summarization function *featureCounts*, which takes as input the SAM or BAM files and assigns them to genomic features. This gives the number of reads mapped per gene, which can then be transformed into RPKM values (Read Per Killobase per Million), normalised and tested for differential expression. 

For the purpose of this course we will use the index previously built and its corresponding GTF file to identify and count the reads mapped to each feature (gene).


```r
# Getting read counts using the index previously built
mycounts<-featureCounts(outputsamfile, annot.ext=file.path(RSUBREAD_INDEX_PATH,"hg19.genes.gtf"), isGTFAnnotationFile=TRUE, isPairedEnd=TRUE)

# Checking your output object
summary(mycounts)
dim(mycounts$counts)
head(mycounts$annotation)
mycounts$targets
mycounts$stat
```



> ## NOTE: {.callout}
>
>  In-built annotations for mm9, mm10 and hg19 are also provided for users' convenience, but we won't be using it for this course.
> **Please do not run the command below**
>
> Getting read counts using the inbuilt hg19 annotation:
>
> mycounts.inbuilt <- featureCounts(myoutputfile, annot.inbuilt="hg19", isGTFAnnotationFile=FALSE, isPairedEnd=TRUE)
> 

For the purpose of this course the read summarisation step has already been performed for all libraries. You will need to load the corresponding RData file to get these read counts. 



```r
MAPPINGDIR <- "/mnt/RNAseqCourse/mapping"
# load the counts previously calculated
load(file.path(MAPPINGDIR,"RawCounts.RData"))
```

```
## Warning in readChar(con, 5L, useBytes = TRUE): cannot open compressed file
## '/mnt/RNAseqCourse/mapping/RawCounts.RData', probable reason 'No such file
## or directory'
```

```
## Error in readChar(con, 5L, useBytes = TRUE): cannot open the connection
```

```r
# check the presence of read counts for the 8 libraries
summary(counts)
```

```
## Error in summary(counts): object 'counts' not found
```



```r
counts$targets
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

You can then print out these counts in a text file for future reference.


```r
# print out counts table for every sample
write.table(counts$counts,file="~/raw_read_counts.txt",sep="\t", quote=F,append=F)
```

```
## Error in is.data.frame(x): object 'counts' not found
```


> ## NOTE: {.callout}
>
> To save disk space on your machine, remember to convert your SAM files to BAM format using samtools and to compress the fastq files. This can be done out of Rstudio, in a shell window:
>
> ` samtools view -Shb myfile.sam -o myfile.bam `
>
> ` gzip myfile.fq `
>
> Visualisation of the BAM files you created can also be done via a genome browser such as IGV (<http://www.broadinstitute.org/igv/>) after sorting and indexing of those files.
>


### QC and stats

Rsubread provides the number of reads mapped to each gene which can then be used for ploting quality control figures and for differential expression analysis.

QC figures of the mapped read counts can be plotted and investigated for potential outlier libraries and to confirm grouping of samples. 

Before plotting QC figures it is useful to get the experiment design. This will allow labeling the data with the sample groups they belong to, or any other parameter of interest.

The experiment design file corresponding to this study has been downloaded from the ArrayExpress webpage and formatted as a tab separated file for this analysis purposes. You can find it in the shared data directory.



```r
# define the experiment design file (tab separated text file is best)
EXPMT_DESIGN_FILE <- file.path(RNAseqDATADIR,'experiment_design.txt')
# read the experiment design file and save it into memory
experiment_design<-read.table(EXPMT_DESIGN_FILE,header=T,sep="\t")
```

```
## Warning in file(file, "rt"): cannot open file
## '/mnt/RNAseqCourse/raw_data/experiment_design.txt': No such file or
## directory
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
#
# set the rownames to the sampleID to allow for ordering
rownames(experiment_design) <- experiment_design$SampleID
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design' not found
```

```r
# order the design following the counts sample order
experiment_design.ord <- experiment_design[colnames(counts$counts),]
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design' not found
```

```r
# look at the design
experiment_design.ord
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
# list the ordered samples for future use
samples <- as.character(experiment_design.ord$SampleID)
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
# create factors for future plotting
group<-factor(experiment_design.ord$tissue)
```

```
## Error in factor(experiment_design.ord$tissue): object 'experiment_design.ord' not found
```

```r
group
```

```
## Error in eval(expr, envir, enclos): object 'group' not found
```

```r
age<-factor(experiment_design.ord$age)
```

```
## Error in factor(experiment_design.ord$age): object 'experiment_design.ord' not found
```

```r
age
```

```
## Error in eval(expr, envir, enclos): object 'age' not found
```


#### Basic QC plots

Density plots of log-intensity distribution of each library can be superposed on a single graph for a better comparison between libraries and for identification of libraries with weird distribution. On the boxplots the density distributions of raw log-intensities are not expected to be identical but still not totally different.

> ## NOTE: {.callout}
>
> The *png* function will create a png file where to save the plots created straight after, and will close this file when *dev.off()* is called.
>
> To see your plots interactively, simply ommit those two lines.



```r
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
```


<img src="fig/RNAseq/Raw_read_counts_per_gene.density.png" alt="boxplotlog10" />


Boxplots of the raw read counts after log10 transformation.


```r
## box plots of raw read counts (log10)
png(file="~/Raw_read_counts_per_gene.boxplot.png")
logcounts <- log(counts$counts,10)
boxplot(logcounts, main="", xlab="", ylab="Raw read counts per gene (log10)",axes=FALSE)
axis(2)
axis(1,at=c(1:length(samples)),labels=colnames(logcounts),las=2,cex.axis=0.8)
dev.off()
```


<img src="fig/RNAseq/Raw_read_counts_per_gene.boxplot.png" alt="boxplotlog10" />


In order to investigate the relationship between samples, hierarchical clustering can be performed using the *heatmap* function from the *stats* package. In this example *heatmap* calculates a matrix of euclidean distances from the mapped read counts for the 100 most highly expressed genes. 



```r
# select data for the 100 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:100]
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
# heatmap with sample name on X-axis
png(file="~/High_expr_genes.heatmap.png")
heatmap(highexprgenes_counts, col=topo.colors(50), margin=c(10,6))
```

```
## Error in heatmap(highexprgenes_counts, col = topo.colors(50), margin = c(10, : object 'highexprgenes_counts' not found
```

```r
dev.off()
```


You will notice that the samples clustering does not follow the original order in the data matrix (alphabetical order "ERR420386" to "ERR420393"). They have been re-ordered according to the similarity of the 100 genes expression profiles.
To understand what biological effect lies under this clustering, one can use the samples annotation for labeling (samples group, age, sex etc).



```r
# heatmap with condition group as labels
colnames(highexprgenes_counts)<- group
```

```
## Error in eval(expr, envir, enclos): object 'group' not found
```

```r
# plot
png(file="~/High_exprs_genes.heatmap.group.png")
heatmap(highexprgenes_counts, col = topo.colors(50), margin=c(10,6))
```

```
## Error in heatmap(highexprgenes_counts, col = topo.colors(50), margin = c(10, : object 'highexprgenes_counts' not found
```

```r
dev.off()
```

<img src="fig/RNAseq/High_expr_genes.heatmap.group.png" alt="heatmap" />


> ## Exercise: Heatmap  {.challenge}
> Produce a heatmap for the 50 most highly expressed genes and annotate the samples with with their age.
>
> * Subset the read counts object for the 50 most highly expressed genes
> * Annotate the samples in the subset with their age (check order with design!)
> * Plot a heatmap with this subset of data, scaling genes and ordering both genes and samples


#### Principal Component Analysis

A Principal Component Analysis (PCA) can also be performed with these data using the *cmdscale* function (from the *stats* package) which performs a classical multidimensional scaling of a data matrix. 


Reads counts need to be transposed before being analysed with the *cmdscale* functions, i.e. genes should be in columns and samples should be in rows. This is the code for transposing and checking the data before further steps:


```r
# select data for the 1000 most highly expressed genes
select = order(rowMeans(counts$counts), decreasing=TRUE)[1:100]
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
# annotate the data with condition group as labels
colnames(highexprgenes_counts)<- group
```

```
## Error in eval(expr, envir, enclos): object 'group' not found
```

```r
# transpose the data to have variables (genes) as columns
data_for_PCA <- t(highexprgenes_counts)
```

```
## Error in t(highexprgenes_counts): object 'highexprgenes_counts' not found
```

```r
dim(data_for_PCA)
```

```
## Error in eval(expr, envir, enclos): object 'data_for_PCA' not found
```


The *cmdscale* function will calculate a matrix of dissimilarities from your transposed data and will also provide information about the proportion of explained variance by calculating Eigen values.


```r
## calculate MDS (matrix of dissimilarities)
mds <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)  
```

```
## Error in as.matrix(x): object 'data_for_PCA' not found
```

```r
# k = the maximum dimension of the space which the data are to be represented in
# eig = indicates whether eigenvalues should be returned
```

The variable mds\$eig provides the Eigen values for the first 8 principal components:


```r
mds$eig
```

```
## Error in eval(expr, envir, enclos): object 'mds' not found
```

Plotting this variable as a percentage will help you determine how many components can explain the variability in your dataset and thus how many dimensions you should be looking at.


```r
# transform the Eigen values into percentage
eig_pc <- mds$eig * 100 / max(mds$eig)
# plot the PCA
png(file="~/PCA_PropExplainedVariance.png")
plot(eig_pc,
     type="h", lwd=15, las=1,
     xlab="Dimensions", 
     ylab="Proportion of explained variance", y.axis=NULL,
     col="darkgrey")
dev.off()
```


<img src="fig/RNAseq/PCA_PropExplainedVariance.png" alt="prop" />


In most cases, the first 2 components explain more than half the variability in the dataset and can be used for plotting. The *cmdscale* function run with default parameters will perform a principal components analysis on the given data matrix and the *plot* function will provide scatter plots for individuals representation. 


```r
## calculate MDS
mds <- cmdscale(dist(data_for_PCA)) # Performs MDS analysis 
```

```
## Error in as.matrix(x): object 'data_for_PCA' not found
```



```r
#Samples representation
png(file="~/PCA_Dim1vsDim2.png")
plot(mds[,1], -mds[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(mds[,1], -mds[,2], rownames(mds), cex=0.8) 
dev.off()
```


<img src="fig/RNAseq/PCA_Dim1vsDim2.group.png" alt="PCA1" />

The PCA plot of the first two components show a clear separation of the Brain and Liver samples across the 1st dimension. Within each sample group we can also notice a split between the 4 samples of each group, which seem to cluster in pair. This observation can be explained by another factor of variability in the data, commonly batch effect or another biological biais such as age or sex.


> ## Exercise: PCA {.challenge}
> Produce a PCA plot from the read counts of the 500 most highly expressed genes and change the labels until you can identify the reason for the split between samples from the same tissue. 
>
> * Get the read counts for the 500 most highly expressed genes
> * Transpose this matrix of read counts
> * Check the number of dimensions explaining the variability in the dataset
> * Run the PCA with an appropriate number of components
> * Annotate the samples with their age \& re-run the PCA \& plot the main components
> * Annotate the samples with other clinical data \& re-run the PCA \& plot the main components until you can separate the samples within each tissue group





### Differential Expression


Before proceeding with differential expression analysis, it is useful to filter out very lowly expressed genes. This will help increasing the statistical power of the analysi while keeping genes of interest. A common way to do this is by filtering out genes having less than 1 count-per-million reads (cpm) in half the samples. 
The "edgeR" library provides the *cpm* function which can be used here.


```r
# load required libraries
library(edgeR)
```

```
## Loading required package: methods
```

```
## Loading required package: limma
```


```r
# get the expression counts from previous alignment step
mycounts <- counts$counts
```

```
## Error in eval(expr, envir, enclos): object 'counts' not found
```

```r
dim(mycounts)
```

```
## Error in eval(expr, envir, enclos): object 'mycounts' not found
```

```r
mycounts[1:5,1:3]
```

```
## Error in eval(expr, envir, enclos): object 'mycounts' not found
```

```r
# filtering
#Keep genes with least 1 count-per-million reads (cpm) in at least 4 samples
isexpr <- rowSums(cpm(mycounts)>1) >= 4
```

```
## Error in cpm(mycounts): object 'mycounts' not found
```

```r
table(isexpr)
```

```
## Error in table(isexpr): object 'isexpr' not found
```

```r
mycounts <- mycounts[isexpr,]
```

```
## Error in eval(expr, envir, enclos): object 'mycounts' not found
```

```r
genes <- rownames(mycounts)
```

```
## Error in rownames(mycounts): object 'mycounts' not found
```

```r
dim(mycounts)
```

```
## Error in eval(expr, envir, enclos): object 'mycounts' not found
```



The *limma* package (since version 3.16.0) offers the *voom* function that will normalise read counts and apply a linear model to the normalised data before computing moderated t-statistics of differential expression.



```r
# load required libraries
library(limma)

# check your samples grouping
experiment_design.ord[colnames(mycounts),]$tissue == group
```

```
## Error in eval(expr, envir, enclos): object 'experiment_design.ord' not found
```

```r
# create design matrix for limma
design <- model.matrix(~0+group)
```

```
## Error in eval(expr, envir, enclos): object 'group' not found
```

```r
# substitute "group" from the design column names
colnames(design)<- gsub("group","",colnames(design))
```

```
## Error in is.data.frame(x): object 'design' not found
```

```r
# check your design matrix
design
```

```
## Error in eval(expr, envir, enclos): object 'design' not found
```

```r
# calculate normalization factors between libraries
nf <- calcNormFactors(mycounts)
```

```
## Error in is(object, "DGEList"): object 'mycounts' not found
```

```r
#
# normalise the read counts with 'voom' function
y <- voom(mycounts,design,lib.size=colSums(mycounts)*nf)
```

```
## Error in is(counts, "DGEList"): object 'mycounts' not found
```

```r
# extract the normalised read counts
counts.voom <- y$E
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
# save normalised expression data into output dir
write.table(counts.voom,file="~/counts.voom.txt",row.names=T,quote=F,sep="\t")
```

```
## Error in is.data.frame(x): object 'counts.voom' not found
```

```r
# fit linear model for each gene given a series of libraries
fit <- lmFit(y,design)
```

```
## Error in is(object, "list"): object 'y' not found
```

```r
# construct the contrast matrix corresponding to specified contrasts of a set of parameters
cont.matrix <- makeContrasts(liver-brain,levels=design)
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
# compute estimated coefficients and standard errors for a given set of contrasts
fit <- contrasts.fit(fit, cont.matrix)
```

```
## Error in NCOL(fit$coefficients): object 'fit' not found
```

```r
# compute moderated t-statistics of differential expression by empirical Bayes moderation of the standard errors
fit <- eBayes(fit)
```

```
## Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, : object 'fit' not found
```

```r
options(digits=3)

# check the output fit
dim(fit)
```

```
## Error in eval(expr, envir, enclos): object 'fit' not found
```

The*topTable* function summarises the output from limma in a table format. Significant DE genes for a particular comparison can be identified by selecting genes with a p-value smaller than a chosen cut-off value and/or a fold change greater than a chosen value in this table. By default the table will be sorted by increasing adjusted p-value, showing the most significant DE genes at the top.



```r
# set adjusted pvalue threshold and log fold change threshold
mypval=0.01
myfc=3

# get the coefficient name for the comparison  of interest
colnames(fit$coefficients)
```

```
## Error in is.data.frame(x): object 'fit' not found
```

```r
mycoef="liver - brain"
# get the output table for the 10 most significant DE genes for this comparison
topTable(fit,coef=mycoef)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
# get the full table ("n = number of genes in the fit")
limma.res <- topTable(fit,coef=mycoef,n=dim(fit)[1])
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
# get significant DE genes only (adjusted p-value < mypval)
limma.res.pval <- topTable(fit,coef=mycoef,n=dim(fit)[1],p.val=mypval)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
dim(limma.res.pval)
```

```
## Error in eval(expr, envir, enclos): object 'limma.res.pval' not found
```

```r
# get significant DE genes with low adjusted p-value high fold change
limma.res.pval.FC <- limma.res.pval[which(abs(limma.res.pval$logFC)>myfc),]
```

```
## Error in eval(expr, envir, enclos): object 'limma.res.pval' not found
```

```r
dim(limma.res.pval.FC)
```

```
## Error in eval(expr, envir, enclos): object 'limma.res.pval.FC' not found
```

```r
# write limma output table for significant genes into a tab delimited file
filename = paste("~/DEgenes_",mycoef,"_pval",mypval,"_logFC",myfc,".txt",sep="")
write.table(limma.res.pval.FC,file=filename,row.names=T,quote=F,sep="\t")
```

```
## Error in is.data.frame(x): object 'limma.res.pval.FC' not found
```

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





### Gene Annotation 

The annotation of EntrezGene IDs from RNAseq data can be done using the BioMart database which contains many species including Human, Mouse, Zebrafish, Chicken and Rat.




```r
# get the Ensembl annotation for human genome
library(biomaRt)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL",host="www.ensembl.org"))

# get entrez gene IDs from limma output table
entrez_genes <- as.character(rownames(limma.res.pval.FC))
```

```
## Error in rownames(limma.res.pval.FC): object 'limma.res.pval.FC' not found
```

```r
length(entrez_genes)
```

```
## Error in eval(expr, envir, enclos): object 'entrez_genes' not found
```

```r
# interrogate the BioMart database to get gene symbol and description for these genes 
detags.IDs <- getBM(
 filters= "entrezgene",
 attributes= c("entrezgene","hgnc_symbol","description"),
 values= entrez_genes,
 mart= mart
)
```

```
## Error in getBM(filters = "entrezgene", attributes = c("entrezgene", "hgnc_symbol", : object 'entrez_genes' not found
```

```r
dim(detags.IDs)
```

```
## Error in eval(expr, envir, enclos): object 'detags.IDs' not found
```

```r
head(detags.IDs)
```

```
## Error in head(detags.IDs): object 'detags.IDs' not found
```



In many cases, several annotations ar available per entrez gene ID. This results in duplicate entries in the output table from *getBM*. The simplest way to deal with this issue is to remove duplicates, although they can also be concatenated in some ways.

Once the annotation has been obtained for all DE genes, this table can be merged with the output table from limma for a complete result and an easier interpretation.



```r
# remove duplicates
detags.IDs.matrix<-detags.IDs[-which(duplicated(detags.IDs$entrezgene)),]
```

```
## Error in eval(expr, envir, enclos): object 'detags.IDs' not found
```

```r
# select genes of interest only
rownames(detags.IDs.matrix)<-detags.IDs.matrix$entrezgene
```

```
## Error in eval(expr, envir, enclos): object 'detags.IDs.matrix' not found
```

```r
entrez_genes.annot <- detags.IDs.matrix[as.character(entrez_genes),]
```

```
## Error in eval(expr, envir, enclos): object 'detags.IDs.matrix' not found
```

```r
# join the two tables
rownames(limma.res.pval.FC) <- limma.res.pval.FC$ID
```

```
## Error in eval(expr, envir, enclos): object 'limma.res.pval.FC' not found
```

```r
limma.res.pval.FC.annot <- cbind(entrez_genes.annot,limma.res.pval.FC)
```

```
## Error in cbind(entrez_genes.annot, limma.res.pval.FC): object 'entrez_genes.annot' not found
```

```r
# check the annotated table
head(limma.res.pval.FC.annot)
```

```
## Error in head(limma.res.pval.FC.annot): object 'limma.res.pval.FC.annot' not found
```



### Gene Set Enrichment

Gene Ontology (GO) enrichment is a method for investigating sets of genes using the Gene Ontology system of classification, in which genes are assigned to a particular set of terms for three major domains: cellular component, biological process and molecular function.

The *GOstats* package can test for both over and under representation of GO terms using the Hypergeometric test. The output of the analysis is typically a ranked list of GO terms, each associated with a p-value.

The Hypergeometric test will require both a list of selected genes (i.e. your DE genes) and a "universe" list (e.g. all genes annotated in the genome you are working with), all represented by their "EntrezGene" ID.


```r
# load the library
library(GOstats)
```

```
## Error in library(GOstats): there is no package called 'GOstats'
```

```r
# Define list of genes of interest (DE genes - EntrezGene IDs)
entrezgeneids <- as.character(rownames(limma.res.pval.FC))
```

```
## Error in rownames(limma.res.pval.FC): object 'limma.res.pval.FC' not found
```

```r
length(entrezgeneids)
```

```
## Error in eval(expr, envir, enclos): object 'entrezgeneids' not found
```

```r
# Define the universe
universeids <- rownames(mycounts)
```

```
## Error in rownames(mycounts): object 'mycounts' not found
```

```r
length(universeids)
```

```
## Error in eval(expr, envir, enclos): object 'universeids' not found
```



Before running the hypergeometric test with the *hyperGTest* function, you need to define the parameters for the test (gene lists, ontology, test direction) as well as the annotation database to be used. The ontology to be tested can be any of the three GO domains: biological process ("BP"), cellular component ("CC") or molecular function ("MF").

In the example below we will test for over-represented biological processes in our list of differentially expressed genes.



```r
# define the p-value cut off for the hypergeometric test
hgCutoff <- 0.05
params <- new("GOHyperGParams",annotation="org.Hs.eg",geneIds=entrezgeneids,universeGeneIds=universeids,ontology="BP",pvalueCutoff=hgCutoff,testDirection="over")
```

```
## Error in getClass(Class, where = topenv(parent.frame())): "GOHyperGParams" is not a defined class
```

```r
#  Run the test
hg <- hyperGTest(params)
```

```
## Error in eval(expr, envir, enclos): could not find function "hyperGTest"
```

```r
# Check results
hg
```

```
## Error in eval(expr, envir, enclos): object 'hg' not found
```



You can get the output table from the test for significant GO terms only by adjusting the pvalues with the *p.adjust* function.


```r
## Get the p-values of the test
hg.pv <- pvalues(hg)
```

```
## Error in eval(expr, envir, enclos): could not find function "pvalues"
```

```r
## Adjust p-values for multiple test (FDR)
hg.pv.fdr <- p.adjust(hg.pv,'fdr')
```

```
## Error in p.adjust(hg.pv, "fdr"): object 'hg.pv' not found
```

```r
## select the GO terms with adjusted p-value less than the cut off
sigGO.ID <- names(hg.pv.fdr[hg.pv.fdr < hgCutoff])
```

```
## Error in eval(expr, envir, enclos): object 'hg.pv.fdr' not found
```

```r
length(sigGO.ID)
```

```
## Error in eval(expr, envir, enclos): object 'sigGO.ID' not found
```

```r
# get table from HyperG test result
df <- summary(hg)
```

```
## Error in summary(hg): object 'hg' not found
```

```r
# keep only significant GO terms in the table
GOannot.table <- df[df[,1] %in% sigGO.ID,]
```

```
## Error in df[, 1]: object of type 'closure' is not subsettable
```

```r
head(GOannot.table)
```

```
## Error in head(GOannot.table): object 'GOannot.table' not found
```


The Gene Ontology enrichment result can be saved in a text file or an html file for future reference.


```r
# Create text report of the significantly over-represented GO terms
write.table(GOannot.table,file="~/GOterms_OverRep_BP.txt",sep="\t",row.names=F)
```

```
## Error in is.data.frame(x): object 'GOannot.table' not found
```

```r
# Create html report of all over-represented GO terms
htmlReport(hg, file="~/GOterms_OverRep_BP.html")
```

```
## Error in eval(expr, envir, enclos): could not find function "htmlReport"
```


> ## Exercise: GOstats {.challenge}
> Identify the GO terms in the Molecular Function domain that are over-represented (pvalue<0.01) in your list of DE genes.
>
> * Get your list of DE genes (Entrez Gene IDs)
> * Set the new parameters for the hypergeometric test
> * Run the test and adjust the pvalues in the output object
> * Identify the significant GO terms at pvalue 0.01



Other softwares can be used to investigate over-represented pathways, such as GeneGO (<https://portal.genego.com/>) and Ingenuity (<http://www.ingenuity.com/products/ipa>). The advantage of these softwares is that they maintain curated and up-to-date extensive databases. They also provide intuitive visualisation and network modelling tools. 


You can save an image of your RNAseq analysis before moving on to the next part of this course.

```r
RDataFile <- "~/RNAseq_DE_analysis_with_R.RData"
save.image(RDataFile)
```



### Record package and version info with `sessionInfo()`

The *sessionInfo* function prints version information about R and any attached packages. It's a good practice to always run this command at the end of your R session and record it for the sake of reproducibility in the future.


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
## [1] biomaRt_2.18.0 edgeR_3.4.2    limma_3.16.7  
## 
## loaded via a namespace (and not attached):
## [1] formatR_1.4    tools_3.2.0    RCurl_1.95-4.1 knitr_1.13    
## [5] stringr_0.6.2  XML_3.98-1.1   evaluate_0.9
```

