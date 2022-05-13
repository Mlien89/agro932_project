#######################################
# Author: Mason Lien                  #
# Description: Explore GBS Data       #
#######################################


#load packages
library(vcfR)
library(ggplot2)
library(adegenet)
library(poppr)

#load variant call file
vcf <- read.vcfR("C:/Users/u942451/OneDrive - University of Nebraska-Lincoln/school/PhD/Genomic Data/NIN18_data_biallelic.recode.vcf", verbose = F)
vcf1 <- read.vcfR("C:/Users/u942451/OneDrive - University of Nebraska-Lincoln/school/PhD/Genomic Data/NIN18_qc_out_vcf.vcf", verbose = F)

#print summary
vcf
vcf1
strwrap(vcf@meta[1:7])

#observe fix region
head(getFIX(vcf))

#observe gt region
vcf@gt[1:6,1:4]

####QC - Quantify Missing data###
#extract the depth (DP)
dp <- extract.gt(vcf, as.numeric = T)

#calculate total average missing data
total_avg_missing <- sum(is.na(dp[,1]))/dim(dp)[1]

#quantify missingness in all samples
GBS_missing <- apply(dp, MARGIN = 2, function(x){ sum(is.na(x)) })
GBS_missing <- GBS_missing/nrow(dp)

#observe missing-ness visually Pre-QC

#heatmap
heatmap.bp(dp[1:1000,], rlabels = F)

#create bar plot of missing data
bp_missing <- as.data.frame(GBS_missing)
bp_missing$genotype <- row.names(bp_missing)

bp_missing$genotype <- reorder(bp_missing$genotype,-bp_missing$GBS_missing)

NIN_2018 <- 
  ggplot(bp_missing,aes(genotype,bp_missing$GBS_missing))+
  geom_bar(stat="identity")+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))

print(NIN_2018 + labs(y = "% Missing", x = "Genotype"))

#quantify missingness with histogram
hist_missing <- apply(dp, MARGIN = 1, function(x){ sum(is.na(x)) }) 
hist_missing <- hist_missing/ncol(vcf@gt[,-1])
hist(hist_missing, col = "#8DD3C7", xlab = "Missingness (%)", main = "2018 NIN GBS Data")

### GBS QC ###

#omit missing samples > 0.7
vcf@gt <- vcf@gt[,c(TRUE, GBS_missing < 0.7)]
#observe change
vcf

#create new depth variable 
dp_omit1 <- extract.gt(vcf, as.numeric=T)

#omit samples with new dp variable
GBS_omit1 <- apply(dp_omit1, MARGIN = 2, function(x){sum(is.na(x))})
GBS_omit1 <- GBS_omit1/nrow(dp_omit1)

#observe changes - notice number of samples has changed along with variants and overall % missing data
vcf

#produce heatmap of omitted samples
heatmap.bp(dp_omit1[1:1000,], rlabels = F)

#omit variants
GBS_omit2 <- apply(dp_omit1, MARGIN = 1, function(x){sum(is.na(x))})
GBS_omit2 <- GBS_omit2/ncol(dp_omit1)
vcf <- vcf[GBS_omit2 < 0.2,]

#observe changes
vcf

#store new depth variable with omitted variants
dp_omit2 <- extract.gt(vcf, as.numeric = T)

#produce heatmap of omitted variants
heatmap.bp(dp_omit2[1:1000,], rlabels = F)

#### Analysis of genome data ###

#make available the vcf file to adegenet packaged using genlight object
x <- vcfR2genlight(vcf)

#print genlight object
x

#A genlight object only supports biallelic, or binary, variants. That is, variants with no more than two alleles. 
#However, variant call format data can include multiple alleles. 
#When we created our genlight object we recieved a warning message indicating that our vcfR object had variants with more than two alleles and that it was being subset to only biallelic variants. 
#This is one of several important differences in how data is handled in VCF data versus genlight objects.
#Another important difference among VCF and genlight data is how the genotypes are stored. 
#In VCF data the alleles are delimited by either a pipe or a forward slash (‘|’, ‘/’ respectively). 
#Because genlight objects only use biallelic loci the genotypes can be recoded as 0, 1 and 2. 
#These correspond to homozygous for the reference or zero allele, heterozygote or homozygous for the first alternate allele. 
#We can validate this by checking a few select genotypes from both the vcfR object and the genlight object.

#vcfr
gt <- extract.gt(vcf, element = "GT")
gt[c(2,6,18),1:3]

#observe genlight matrix
t(as.matrix(x))[c(1,5,17),1:3]

#assign names to x
pop(x) <- x$ind.names
popNames(x)

#change ploidy level to 2
ploidy(x) <- 2

#Principle Component Analysis
x.pca <- glPca(x, nf = 3)
barplot(100*x.pca$eig/sum(x.pca$eig), col = heat.colors(50), main = "PCA Eigenvalues")
title(ylab = "Percent of variance/nexplained", line = 2)
title(xlab = "Eigenvalues", line = 1)


#plot pca results
x.pca.scores <- as.data.frame(x.pca$scores)
x.pca.scores$pop <- pop(x)

set.seed(9)
p <- ggplot(x.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p

#create distance matrix
x.dist <- dist(x)

#create distance tree
x.dist <- poppr::bitwise.dist(x)

#subsetting a vcfR object to 200 random variants
alphabet <- c("A", "B", "C", "D")

vcf[c(1:10),]

subset.1 <- sample(size = 200, x = c(1:nrow(vcf)))
subset.2 <- sample(size = 200, x = c(1:nrow(vcf)))

identical(subset.1, subset.2)

vcf.sub1 <- vcf[subset.1,]
vcf.sub2 <- vcf[subset.2,]

#Creating a list object to save our subsets in.
x.variant.subset <- vector(mode = "list", length = 50)

# Using a for loop to generate 50 subsets of 200 random variants from the rubi.VCF vcfR object.
for (i in 1:50){
  x.variant.subset[[i]] <- vcf[sample(size = 200, x= c(1:nrow(vcf)))]
}

length(x.variant.subset)

head(x.variant.subset, n=2)


#create GenLight object
x.gl.subset <- lapply(x.variant.subset, function(x) suppressWarnings(vcfR2genlight(x)))
for (i in 1:length(x.gl.subset)){
  ploidy(x.gl.subset[[i]]) <- 2
}

# Creating a simple UPGMA tree per object
library(phangorn)
x.trees <- lapply(x.gl.subset, function (x) upgma(bitwise.dist(x)))
class(x.trees) <- "multiPhylo"

# Overlapping the trees
densiTree(x.trees, consensus = tree, scaleX = T, show.tip.label = F, alpha = 0.1)
title(xlab = "Proportion of variants different")


#create tree
tree <- aboot(x, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)
