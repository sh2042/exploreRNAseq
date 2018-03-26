####################################################################################
# title: rnaseq_analyis_deseq
# author: shunjan
# purpose: Script to explore and assess the quality of RNAseq data
####################################################################################
library (magrittr)
library(dplyr)
library(ggplot2)

# read in the data matrix
gene_count_matrix <- read.csv("gene_count_matrix2.csv", header=TRUE, sep=",")
phenoData <- read.csv("TOFall.csv", header=TRUE, sep=",")
phenoData <- phenoData %>% filter(Gender == "M")
samples <- as.vector(phenoData[,1])

samples <- paste0("X", sep="", samples)
samples <- gsub("-", "_", samples)

#deleted VSD sample X002440_2  - extreme outlier
row.names(gene_count_matrix) <- gene_count_matrix$gene_id

gene_count_matrix <- gene_count_matrix %>% select(samples)

rawreadcounts <- gene_count_matrix


# histogram  #########################################################################


names <- c(colnames(rawreadcounts))
x <- barplot(colSums(as.matrix(rawreadcounts))/1000000, main="Total number of reads per sample (million)", names.arg=names,las=2)

#want 5 million to 10 million reads



# transform the data ################################################################
# this is not part of differential expression
# just looking at distribution of the data

pseudoCount = log2(rawreadcounts + 1)

x <- barplot(colSums(as.matrix(pseudoCount)), main="Reads per sample log2+1", names.arg=names,las=2)

ggplot(rawreadcounts, aes(x = X004321_1)) + ylab(expression(log[2](count + 1))) +
  geom_histogram(colour = "white", fill = "#525252", binwidth = 0.6)


#boxplots to view between sample distribution
install.packages("reshape")
library(reshape)

df = melt(pseudoCount, variable_name = "Samples")
df = data.frame(df, Condition = "sample")
ggplot(df, aes(x = Samples, y = value, fill = Condition)) + geom_boxplot() + xlab("") +
  ylab(expression(log[2](count + 1))) + scale_fill_manual(values = c("#619CFF", "#F564E3"))

# density plots

ggplot(df, aes(x = value, colour = Samples, fill = Samples)) + ylim(c(0, 0.25)) +
  geom_density(alpha = 0.2, size = 1.25) + facet_wrap(~ Condition) +
  theme(legend.position = "top") + xlab(expression(log[2](count + 1)))


#heatmap  not very informative
#install.packages("mixOmics", dependencies = TRUE)
library(mixOmics)  ###### not working
#install.packages("RColorBrewer")
library(RColorBrewer)


mat.dist = pseudoCount
colnames(mat.dist) = paste(colnames(mat.dist), sep = " : ")
mat.dist = as.matrix(dist(t(mat.dist)))
mat.dist = mat.dist/max(mat.dist)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(16)
cim(mat.dist, symkey = FALSE, margins = c(9, 9))



##############################################################################################
# DEseq analysis - differential gene expression

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(DESeq2)

# build the dds

cts <- as.matrix(rawreadcounts)
# include phenodata

coldata <- phenoData
row.names(coldata) <- coldata$Sample.ID
coldata <- coldata[,-1]

# check is coldDAta and cts samples are in order
rownames(coldata) <- paste0("X", sep="", rownames(coldata))
rownames(coldata) <- gsub("-", "_", rownames(coldata))

all(rownames(coldata) == colnames(cts))

cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 1)
dds

########## pre-filtering ######################
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


###### differential expression  ########
dds <- DESeq(dds)
res <- results(dds)
res

res <- results(dds, addMLE=TRUE)
plotMA(res)
plotMA(res, MLE=TRUE)

resultsNames(dds)
resLFC <- lfcshrink(dds, coef=2)
resLFC

#order results by smallest p-value
resOrdered <- res[order(res$pvalue),]


#how many adjusted p values less than one
sum(res$padj < 0.1, na.rm=TRUE)

res05 <- results(dds, alpha=0.05)
summary(res05)
sum(res05$padj < 0.05, na.rm=TRUE)

summary(res)

# examin the read counts across genes
plotCounts(dds, gene="ENSG00000182197", intgroup="Diagnosis")


d <- plotCounts(dds, gene="ENSG00000182197", intgroup="Tissue.Type",
                returnData=TRUE)
library("ggplot2")
ggplot(d, aes(x=Tissue.Type, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))


############    transform data        ################################################

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)  #used rld for analysis


# this gives log2(n + 1)
ntd <- normTransform(dds)
install.packages("vsn")
library(vsn)
plot(assay(rld))


############    visualizing   ################################################

################   heatmap           #########################################
#install.packages("pheatmap")
library(pheatmap)

# heatmap for the count matrix

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Tissue.Type","Diagnosis", "Age")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(rld)))

#heatmap for the sample distance matrix
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, rld$type, sep="-")
#colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "YlOrRd")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         annotation_col=df,
         col=colors)

###### custom heatmap ##########

# explore the count matrix

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Tissue.Type","Diagnosis")])
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)

########## pca plot ##########################

plotPCA(rld, intgroup=c("Age", "Diagnosis"))
# get a plot with just the sample names


which(sampleDistMatrix == max(sampleDistMatrix), arr.ind = TRUE)
apply(sampleDistMatrix,2,max)
max(colSums(sampleDistMatrix))
sampleDistMatrix[,"X07520_RV"]



###########################################################################
## top variant genes #####################################################
#############################################################################

library("genefilter")
topVarGenes <- head(order(-rowVars(assay(rld))),20)
#topVarGenes[ENSG00000182197]
mat <- assay(rld)[ topVarGenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,c("Age", "Diagnosis")])
pheatmap(mat)

#######################################################################
## look at one gene #########################################################
# EXT1, DHCR7, F7, SLC25A1, MYLK, AP3B1, CITED2, PLA2G6


assay(rld)["ENSG00000141510",]
data <- assay(rld)[c("ENSG00000182197", "ENSG00000172893", "ENSG00000057593", "ENSG00000100075", "ENSG00000065534", "ENSG00000132842", "ENSG00000164442", "ENSG00000184381"),]


#ENSG00000182197 ext1
names <- colnames(assay(rld))
barplot(data[8,], main="EXT1", names.arg=names, ylab="log fold change", las=2)
barplot(data,legend.text = rownames(data),beside=T, ylab='Logfold Change', las=2)
