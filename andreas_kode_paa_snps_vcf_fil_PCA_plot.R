
# Set working directory 
setwd("~/Documents/Bioinformatics/Rclavata")

wd00 <- "/home/hal9000/Documents/shrfldubuntu18/soemrokke"
setwd(wd00)

wd03d <- "03d_gstacks"
# define input file directory
inf01 <- "03d_gstacks"
# paste together path and input files
pthinf01 <- paste0(wd00,"/",wd03d,"/","populations.haps.vcf")
pthinf02 <- paste0(wd00,"/",wd03d,"/","populations.snps.vcf")


# Loading packages
library(readxl)
library(vcfR)
library(poppr)
library(RColorBrewer)
library(ggplot2)

# Loading VCF files
clavata.snps <- read.vcfR("populations.snps.vcf")
clavata.snps <- read.vcfR(pthinf02)

# Loading population data 
pop.data <- read.table("clavata_popmap.txt", sep = "\t", header = F)
pop.data <- read.table("part07B_popmap.txt", sep = "\t", header = F)

# Naming columns in population data 
colnames(pop.data) <- c("sampleID","location")

# Checking that pop.data and the loaded VCF files contain the same sample names. 
all(colnames(clavata.snps@gt)[-1] == pop.data$sampleID)
## Should return TRUE, if they match. 

# Converting the data set to a genlight object
clavata.snps.gl <- vcfR2genlight(clavata.snps)

# The ploidy of Raja clavata is diploid, so we place a 2 in the slot of the genlight object.
ploidy(clavata.snps.gl) <- 2

# We insert our population data into the genlight object. 
pop(clavata.snps.gl) <- pop.data$location

# Making PCA ---- 
## Creating a table with eigenvalues
clavata.pca <- glPca(clavata.snps.gl, nf = 3)

# Creating the barplot
barplot(100*clavata.pca$eig/sum(clavata.pca$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)

# Extracting the PCA scores from the principal components
clavata.pca.scores <- as.data.frame(clavata.pca$scores)
clavata.pca.scores$pop <- pop(clavata.snps.gl)

#creating the PCA
PCA <- ggplot(clavata.pca.scores, aes(x=PC1, y=PC2, colour=pop)) +
  geom_point(size=2)+
  theme_bw()

set.seed(9)
PCA1 <- ggplot(clavata.pca.scores, aes(x=PC1, y=PC3, colour=pop)) +
  geom_point(size=2)+
  theme_bw()


