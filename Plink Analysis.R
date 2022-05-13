setwd("C:/Users/u942451/OneDrive - University of Nebraska-Lincoln/school/PhD/Genomic Data/")


#recode vcf file into plink for down stream analysis
system("plink.exe --const-fid 0 --vcf NIN18_data.recode.vcf --recode --nonfounders --allow-no-sex --allow-extra-chr --out NIN18")

#.bed formatted plink file
system("plink.exe --file NIN18 --make-bed --nonfounders --allow-no-sex --allow-extra-chr --out NIN18_bed")

#QC options

#1. missingness per SNP?: --geno
#2. missingness per individual: --mind
#3. minor allele frequency: --maf

system("plink.exe --bfile NIN18_bed --geno 0.2 --mind 0.7 --maf 0.05 --make-bed --allow-extra-chr --out NIN18_qc_out")

#make new vcf file
system("plink.exe --bfile NIN18_qc_out --recode vcf-iid --allow-extra-chr --out NIN18_qc_out_vcf")

# PCA analysis: a model-free method ------------------------------------------------------------------

#identify prune sites

system("plink.exe --vcf NIN18_qc_out_vcf.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 2 0.2 --out NIN18")


#pca analysis
system("plink.exe --vcf NIN18_qc_out_vcf.vcf --double-id --allow-extra-chr --set-missing-var-ids @:# --extract NIN18.prune.in --pca --make-bed --out NIN18_prune_pca")


#pca plot
library(tidyverse)

plinkPCA <- read_table2("NIN18_prune_pca.eigenvec", col_names = F)
plinkPCA <- plinkPCA[,c(-1,-2)]
EigenValue <- scan("NIN18_prune_pca.eigenval")

#set column names
names(plinkPCA)[1:ncol(plinkPCA)] <- paste0("PC", 1:(ncol(plinkPCA)))

#percentage variance explained
pve <- data.frame(PC = 1:20, pve = EigenValue/sum(EigenValue)*100)

#make the plot
P <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
P + ylab("Percentage variance explained") + theme_light()

#plot pca
ggplot(plinkPCA, aes(PC1, PC2)) + geom_point(size = 3) + coord_equal() + theme_light() + coord_equal() + 
  theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

