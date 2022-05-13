#######################################
# Author: Mason Lien                  #
# Description: Explore GBS Data       #
#######################################

#load packages
library(vcfR)
library(ggplot2)
library(adegenet)
library(poppr)
library(rrBLUP)

setwd("/Users/masonlien/OneDrive - University of Nebraska-Lincoln/school/PhD/")


#load variant call file
vcf <- read.vcfR("Genomic Data/NIN18_data_biallelic.recode.vcf", verbose = F)
vcf1 <- read.vcfR("Genomic Data/NIN18_qc_out_vcf.vcf", verbose = F)

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

#convert to data.frame for further QC
geno <- as.data.frame(x)

#sum NA values
sum(is.na(geno))

#missing rate
missing <- apply(geno,2, function(x){sum(is.na(x))/length(x)})
#calculate maf
maf <- apply(geno,2, function(x){
  frq <- mean(x, na.rm =T)/2 # 1 allele
  return(ifelse(frq > 0.5, 1-frq, frq))
})

#plot the results
hist(missing, breaks=100, col = "blue", xlab = "SNP Missing rate")
hist(maf, breaks=100, col="blue", xlab = "Minor Allele Freq")

#remove snps with high missing rate and low MAF
idx1 <- which(missing > 0.2)
idx2 <- which(maf < 0.05)
idx <- unique(c(idx1, idx2))

geno2 <- geno[, idx]
dim(geno2)

#import the phenotype data
phenotype <- read.csv("/Volumes/ML - External/2017-2018 season/NIN/NIN 18 BLUPS_GP.csv", header = T)
pheno <- phenotype[phenotype$Genotype %in% row.names(geno2),]
pheno1 <- pheno[,c(1:2)]
#sort order of pheno by genotype file
pheno <- pheno[match(row.names(geno2), pheno$Genotype), ]

#remove any genoyptes not part of the phenotype file
geno2 <- geno2[row.names(geno2) %in% pheno$Genotype,]
geno2 <- geno2 -1
geno2 <- cbind(pheno1,geno2)


#impute missing SNP_Markers with A.mat function
SNP_Markers <- geno2[,3:ncol(geno2)]
impute = (A.mat(SNP_Markers,max.missing=0.5,impute.method="mean",return.imputed=T))
SNP_Markers_1_impute <- impute$imputed
SNP_Markers_1_impute[1:5,1:10]


#####select phenotypic data#########
Pheno <- pheno %>% select(NDVI)

#training and testing sets
training_entries <- as.matrix(sample(1:29, 23))
testing_entries <- setdiff(1:29,training_entries)

Pheno_training_data <- as.matrix(Pheno[training_entries,])
SNP_training_data <- as.matrix(SNP_Markers_1_impute[training_entries,])

Pheno_testing_data <- as.matrix(Pheno[testing_entries,])
SNP_testing_data <- as.matrix(SNP_Markers_1_impute[testing_entries,], K=NULL)

trained_model <- mixed.solve(y=Pheno_training_data, Z = SNP_training_data)

marker_effects <- as.matrix(trained_model$u)
BLUE <- as.vector(trained_model$beta)

predicted_train <- as.matrix(SNP_training_data) %*% marker_effects
predicted_test <- as.matrix(SNP_testing_data) %*% marker_effects
predicted_train_result <- as.vector((predicted_train[,1])+BLUE)
predicted_test_result <- as.vector((predicted_test[,1])+BLUE)

summary(as.vector(predicted_train_result))
summary(predicted_test_result)

cor(as.vector(Pheno_testing_data), predicted_test_result, use="complete")
cor(as.vector(Pheno_training_data), predicted_train_result, use = "complete")

plot(Pheno_testing_data,predicted_test_result)
plot(Pheno_training_data,predicted_train_result)
