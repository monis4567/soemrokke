#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:

# GBS Population genetics on thorny skates
# See these webpages for inspiration:
#https://cran.r-project.org/web/packages/vcfR/vignettes/intro_to_vcfR.html
#https://knausb.github.io/vcfR_documentation/visualization_1.html
#https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
library(readxl)
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
# define working directory
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/soemrokke"
setwd(wd00)
wd03d <- "03d_gstacks"
# define input file directory
inf01 <- "03d_gstacks"
# paste together path and input files
pthinf01 <- paste0(wd00,"/",wd03d,"/","populations.haps.vcf")
pthinf02 <- paste0(wd00,"/",wd03d,"/","populations.snps.vcf")


library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/soemrokke"
setwd(wd00)

wd03d <- "03d_gstacks"
# define input file directory
inf01 <- "03d_gstacks"
# paste together path and input files
pthinf01 <- paste0(wd00,"/",wd03d,"/","populations.haps.vcf")
pthinf02 <- paste0(wd00,"/",wd03d,"/","populations.snps.vcf")
haps.VCF <- read.vcfR(pthinf01)
snps.VCF <- read.vcfR(pthinf02)

pop.data <- read.table("part07B_popmap.txt", sep = "\t", header = TRUE)
colnames(pop.data) <- c("AccessID","location")
#one row is missing , make it a vector
P08899mss <- c("P08899","NorthSea_Thornbackskate")
# bind this as an extra row to the data frame
pop.data <- rbind(pop.data,P08899mss)
pop.data$location[pop.data$AccessID=="P08912"] <- "NorthSea_Thornbackskate"
# use gsub to swap first and second around
pop.data$location <- gsub("(.*)_(.*)","\\2_\\1",pop.data$location)

#check for missing samples
missingsmplh <- setdiff(colnames(haps.VCF@gt)[-1], pop.data$AccessID )
missingsmpls <- setdiff(colnames(snps.VCF@gt)[-1], pop.data$AccessID )
# check number of elements, this should match
length(colnames(haps.VCF@gt)[-1])
length(pop.data$AccessID)
#order the pop.data data frame
pop.data <- pop.data[order(pop.data$AccessID),]
pop.data$AccessID <- as.character(pop.data$AccessID)
# use match function to get the same order of AccessID in the pop.data
# data frame as in the haps.VCF object, otherwise the 'all' comparison 
# will fail
pop.data$AccessID <- pop.data$AccessID[match(colnames(haps.VCF@gt)[-1],pop.data$AccessID)]
pop.data$location <- pop.data$location[match(colnames(haps.VCF@gt)[-1],pop.data$AccessID)]
# and for the snps
pop.data.snp <- pop.data
# mach again
pop.data.snp$AccessID <- pop.data$AccessID[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
pop.data.snp$location <- pop.data$location[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
# now check the order is the same -  this should return TRUE for all
identical(colnames(haps.VCF@gt)[-1], pop.data$AccessID) 
all(colnames(haps.VCF@gt)[-1] == pop.data$AccessID)
#Converting the dataset to a genlight object
gl.haps <- vcfR2genlight(haps.VCF)
# A warning is shown while transforming the object, telling us that there 
# are loci with more than two alleles. This is a diploid
# organism, so we will specify a ploidy of two. 
ploidy(gl.haps) <- 2
# Add predetermined populations. 
pop(gl.haps) <- pop.data$location

#Population genetic analyses for GBS data
#Distance matrices
#create a pairwise genetic distance matrix for individuals or populations 
#To summarize, we can create a distance matrix from a genlight object using dist():
#x.dist <- dist(x)
gl.haps.dist <- poppr::bitwise.dist(gl.haps, euclidean = T)
#make it a matrix to make it a data frame
df_gl.haps.dist <- as.data.frame(as.matrix(gl.haps.dist))
#_______________________________________________________________________________
# start - see the matrix coloured 
#_______________________________________________________________________________
library(tidyverse)
library(ggplot2)
# use tidyverse to re arrange dataframe
df_ghd2 <- df_gl.haps.dist %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#plot data frame as heat map
plt_ghd2 <- ggplot(df_ghd2, aes(x = rowname,
                          y = colname, 
                          fill = value)) +
  geom_tile()

#_________________________________________________________________
# end - see the matrix coloured 
#_________________________________________________________________
#

#Distance tree
#build a genetic distance tree 
tree.gl.haps <- aboot(gl.haps, tree = "upgma", 
                      distance = bitwise.dist, sample = 100, 
                      showtree = F, cutoff = 50, quiet = T)
#  color the tips of the tree based on the population of origin of the samples
cols <- brewer.pal(n = nPop(gl.haps), name = "Dark2")
#_______________________________________________________________________________
# cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
#                  "yellow","white")
# colfunc <- colorRampPalette(cbbPalette2)
# #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# #
# cols <- colfunc(nPop(gl.snps))
#_______________________________________________________________________________


plot.phylo(tree.gl.haps, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.haps)])
nodelabels(tree.gl.haps$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
# get unique location names
locat_gl_haps <- unique(pop(gl.haps))
# add a legend
legend('topleft', 
       legend = c(locat_gl_haps), fill = cols, border = FALSE, bty = "n", cex = 1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#
# Make 
# Minimum spanning networks
library(igraph)

haps.dist <- bitwise.dist(gl.haps)
haps.msn <- poppr.msn(gl.haps, haps.dist, showplot = T, 
                      include.ties = T)

node.size <- rep(2, times = nInd(gl.haps))
names(node.size) <- indNames(gl.haps)
vertex.attributes(haps.msn$graph)$size <- node.size


cols <- brewer.pal(n = nPop(gl.haps), name = "Dark2")
#_______________________________________________________________________________
# cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
#                  "yellow","white")
# colfunc <- colorRampPalette(cbbPalette2)
# #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# #
# cols <- colfunc(nPop(gl.haps))
#_______________________________________________________________________________


set.seed(9)
plot_poppr_msn(gl.haps, haps.msn , 
               palette = brewer.pal(n = nPop(gl.haps), 
                                    name = "Dark2"), gadj = 70)
#Principal components analysis
#A principal components analysis (PCA) converts the observed SNP data
#into a set of values of linearly uncorrelated variables called principal 
#components that summarize the variation between samples. 
#We can perform a PCA on our genlight object by using the glPCA function.
# see barplot first
haps.pca <- glPca(gl.haps, nf = 3)
barplot(100*haps.pca$eig/sum(haps.pca$eig),
        col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
# second make PCA plot
haps.pca.scores <- as.data.frame(haps.pca$scores)
haps.pca.scores$pop <- pop(gl.haps)

library(ggplot2)
set.seed(9)
p <- ggplot(haps.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p01 <- p
# DAPC
pnw.dapc <- dapc(gl.haps, n.pca = 3, n.da = 2)
# check if DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, 
        clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)
pnw.dapc <- dapc(gl.haps, n.pca = 3, n.da = 2)
# compoplot
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'top')
# separate the samples by population
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.haps)
dapc.results$indNames <- rownames(dapc.results)
# transform the data frame usingpivot_longer from the package tidyr
# to make it match input format in ggplot
library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
# see header
head(dapc.results, n = 6)
# Rename the columns into more familiar terms:
colnames(dapc.results) <- c("Original_Pop",
                            "Sample",
                            "Assigned_Pop",
                            "Posterior_membership_probability")
# 
# plot the dapc.results data frame reorganized using pivot_longer, 
# using the samples on the X-axis and membership probabilities on the Y-axis. 
# The fill color will indicate the original population assignments. 
# Each facet represents the original population assignment for each sample:
#   
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, 
                              fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p02 <- p
library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- p01 + labs(title = "a")#,
p02t <- p02 + labs(title = "b")#,


pA <-  p01t +
  p02t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig01_pca_soemrokke_haps.png")
figname02 <- paste(wd00,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}

#
#________________________________________________________________________________
#________________________________________________________________________________
#________________________________________________________________________________

haps.VCF <- read.vcfR(pthinf01)
snps.VCF <- read.vcfR(pthinf02)


#check for missing samples
missingsmplh <- setdiff(colnames(haps.VCF@gt)[-1], pop.data$AccessID )
missingsmpls <- setdiff(colnames(snps.VCF@gt)[-1], pop.data$AccessID )
# check number of elements, this should match
length(colnames(snps.VCF@gt)[-1])
length(pop.data$AccessID)
#order the pop.data data frame
pop.data <- pop.data[order(pop.data$AccessID),]
pop.data$AccessID <- as.character(pop.data$AccessID)
# use match function to get the same order of AccessID in the pop.data
# data frame as in the snps.VCF object, otherwise the 'all' comparison 
# will fail
pop.data$AccessID <- pop.data$AccessID[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
pop.data$location <- pop.data$location[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
# and for the snps
pop.data.snp <- pop.data
# mach again
pop.data.snp$AccessID <- pop.data$AccessID[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
pop.data.snp$location <- pop.data$location[match(colnames(snps.VCF@gt)[-1],pop.data$AccessID)]
# now check the order is the same -  this should return TRUE for all
identical(colnames(snps.VCF@gt)[-1], pop.data$AccessID) 
all(colnames(snps.VCF@gt)[-1] == pop.data$AccessID)
#Converting the dataset to a genlight object
gl.snps <- vcfR2genlight(snps.VCF)
# A warning is shown while transforming the object, telling us that there 
# are loci with more than two alleles. This is a diploid
# organism, so we will specify a ploidy of two. 
ploidy(gl.snps) <- 2
# Add predetermined populations. 
pop(gl.snps) <- pop.data$location

#Population genetic analyses for GBS data
#Distance matrices
#create a pairwise genetic distance matrix for individuals or populations 
#To summarize, we can create a distance matrix from a genlight object using dist():
#x.dist <- dist(x)
gl.snps.dist <- poppr::bitwise.dist(gl.snps, euclidean = T)
#make it a matrix to make it a data frame
df_gl.snps.dist <- as.data.frame(as.matrix(gl.snps.dist))
#_______________________________________________________________________________
# start - see the matrix coloured 
#_______________________________________________________________________________
library(tidyverse)
library(ggplot2)
# use tidyverse to re arrange dataframe
df_ghd2 <- df_gl.snps.dist %>%
  rownames_to_column() %>%
  gather(colname, value, -rowname)
#plot data frame as heat map
plt_ghd2 <- ggplot(df_ghd2, aes(x = rowname,
                                y = colname, 
                                fill = value)) +
  geom_tile()

#_________________________________________________________________
# end - see the matrix coloured 
#_________________________________________________________________
#

#Distance tree
#build a genetic distance tree 
tree.gl.snps <- aboot(gl.snps, tree = "upgma", 
                      distance = bitwise.dist, sample = 100, 
                      showtree = F, cutoff = 50, quiet = T)
#  color the tips of the tree based on the population of origin of the samples
cols <- brewer.pal(n = nPop(gl.snps), name = "Dark2")
#_______________________________________________________________________________
# cbbPalette2 <- c("black","purple","blue","green","yellowgreen",
#                  "yellow","white")
# colfunc <- colorRampPalette(cbbPalette2)
# #https://stackoverflow.com/questions/13353213/gradient-of-n-colors-ranging-from-color-1-and-color-2
# #
# cols <- colfunc(nPop(gl.snps))
#_______________________________________________________________________________
plot.phylo(tree.gl.snps, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(gl.snps)])
nodelabels(tree.gl.snps$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
#legend(35,10,c("CA","OR","WA"),cols, border = FALSE, bty = "n")
# get unique location names
locat_gl_snps <- unique(pop(gl.snps))
# add a legend
legend('topleft', 
       legend = c(locat_gl_snps), fill = cols, border = FALSE, bty = "n", cex = 1)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")
#
# Make 
# Minimum spanning networks
library(igraph)

snps.dist <- bitwise.dist(gl.snps)
snps.msn <- poppr.msn(gl.snps, snps.dist, showplot = T, 
                      include.ties = T)

node.size <- rep(2, times = nInd(gl.snps))
names(node.size) <- indNames(gl.snps)
vertex.attributes(snps.msn$graph)$size <- node.size

set.seed(9)
plot_poppr_msn(gl.snps, snps.msn , 
               palette = brewer.pal(n = nPop(gl.snps), 
                                    name = "Dark2"), gadj = 70)

#Principal components analysis
#A principal components analysis (PCA) converts the observed SNP data
#into a set of values of linearly uncorrelated variables called principal 
#components that summarize the variation between samples. 
#We can perform a PCA on our genlight object by using the glPCA function.
# see barplot first
snps.pca <- glPca(gl.snps, nf = 3)
barplot(100*snps.pca$eig/sum(snps.pca$eig),
        col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
# second make PCA plot
snps.pca.scores <- as.data.frame(snps.pca$scores)
snps.pca.scores$pop <- pop(gl.snps)

library(ggplot2)
set.seed(9)
p <- ggplot(snps.pca.scores, aes(x=PC1, y=PC2, colour=pop)) 
p <- p + geom_point(size=2)
p <- p + stat_ellipse(level = 0.95, size = 1)
p <- p + scale_color_manual(values = cols) 
p <- p + geom_hline(yintercept = 0) 
p <- p + geom_vline(xintercept = 0) 
p <- p + theme_bw()
p01 <- p
# DAPC
pnw.dapc <- dapc(gl.snps, n.pca = 3, n.da = 2)
# check if DAPC is similar to the PCA we can plot the data in a scatter plot.
scatter(pnw.dapc, col = cols, cex = 2, legend = TRUE, 
        clabel = F, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75)
pnw.dapc <- dapc(gl.snps, n.pca = 3, n.da = 2)
# compoplot
#compoplot(pnw.dapc,col = function(x) cols, posi = 'top')
compoplot(pnw.dapc,col = cols, posi = 'top')
# separate the samples by population
dapc.results <- as.data.frame(pnw.dapc$posterior)
dapc.results$pop <- pop(gl.snps)
dapc.results$indNames <- rownames(dapc.results)
# transform the data frame usingpivot_longer from the package tidyr
# to make it match input format in ggplot
library(reshape2)
# dapc.results <- melt(dapc.results)
library(tidyr)
dapc.results <- pivot_longer(dapc.results, -c(pop, indNames))
# see header
head(dapc.results, n = 6)
# Rename the columns into more familiar terms:
colnames(dapc.results) <- c("Original_Pop",
                            "Sample",
                            "Assigned_Pop",
                            "Posterior_membership_probability")
# 
# plot the dapc.results data frame reorganized using pivot_longer, 
# using the samples on the X-axis and membership probabilities on the Y-axis. 
# The fill color will indicate the original population assignments. 
# Each facet represents the original population assignment for each sample:
#   
p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, 
                              fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p02 <- p
library(patchwork)
# see this example: https://www.datanovia.com/en/blog/ggplot-title-subtitle-and-caption/
#caption = "Data source: ToothGrowth")
p01t <- p01 + labs(title = "a")#,
p02t <- p02 + labs(title = "b")#,


pA <-  p01t +
  p02t +
  
  plot_layout(nrow=2,byrow=T) + #xlab(xlabel) +
  plot_layout(guides = "collect") +
  plot_annotation(caption=pthinf01) #& theme(legend.position = "bottom")
#p
bSaveFigures=T
#make filename to save plot to
figname01 <- paste0("Fig02_pca_soemrokke_snps.png")
figname02 <- paste(wd00,"/",figname01,sep="")
if(bSaveFigures==T){
  ggsave(pA,file=figname02,width=210,height=297,
         units="mm",dpi=300)
}


#

#
