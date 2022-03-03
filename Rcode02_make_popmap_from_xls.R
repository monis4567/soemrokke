#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# GBS Population genetics on thorny skates
library(readxl)
# define working directory
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/soemrokke"
setwd(wd00)
# define input file
inf01 <- "soemrokke_pop_map_table01.xls"
# paste together path and input flie
pthinf01 <- paste0(wd00,"/",inf01)
# read in excel file as tibble, get specific sheet and skip 20 rows
tibl_inx01 <- readxl::read_xls(pthinf01)
# drop first row
#tibl_inx01 <- tibl_inx01[-1,]
#substitute spaces underscores
#colnames(tibl_inx01) <- gsub(" ","_",colnames(tibl_inx01) )
# replace spaces in names
tibl_inx01$Species <- gsub(" ","_",tibl_inx01$Species)
tibl_inx01$Location <- gsub(" ","_",tibl_inx01$Location)
tibl_inx01$Species <- gsub("\\?","unknown",tibl_inx01$Species)
#unique(tibl_inx01$Species)
#unique(tibl_inx01$Location)
df_p01 <- as.data.frame(tibl_inx01)
# take only column number 1 and 3
df_p02 <- df_p01[,c(1,3)]
# Find the unique elements
unqS <- unique(df_p02$sampleNo)
# find the number of unique elements 
#length(unqS)
#Remove duplicated samples
df_p02 <- df_p02[!duplicated(df_p02[1]),]
# Remove samples that never ended up being sequenced
# sample "P08867" was never sequenced
df_p02 <- df_p02[!df_p02[1]=="P08867",]
# I decided to try and remove sample 'P08953' as I kept getting the error
# 
# Processing sample /groups/hologenomics/phq599/data/soemrokke/03_stacks/03d_gstacks/P08953 [27 of 94]
# Parsing /groups/hologenomics/phq599/data/soemrokke/03_stacks/03d_gstacks/P08953.tags.tsv.gz
# Parsing /groups/hologenomics/phq599/data/soemrokke/03_stacks/03d_gstacks/P08953.snps.tsv.gz
# Parsing /groups/hologenomics/phq599/data/soemrokke/03_stacks/03d_gstacks/P08953.alleles.tsv.gz
# srun: error: node925: task 0: Segmentation fault
# srun: Terminating job step 30484789.0
#
# When I tried running it remotely on the HPC server, so sample "P08953" is excluded here 
# and I got the same error for 'SR23'
# and I got the same error for 'P2397638'
df_p02 <- df_p02[!df_p02[1]=="P08953",]
df_p02 <- df_p02[!df_p02[1]=="SR23",]
df_p02 <- df_p02[!df_p02[1]=="P2397638",]

# with the aim getting an index barcode list as described under section 
# '4.1.2 Specifying the barcodes'
#https://catchenlab.life.illinois.edu/stacks/manual/#prun
outfl1 = "part06C_popmap.txt"
# paste together path and input flie
pthoutf01 <- paste0(wd00,"/",outfl1)
# use tab as separator
write.table(df_p02, file=pthoutf01, sep="\t",
            row.names = F, # do not use row names
            col.names = F, # do not use columns names
            quote = F) # do not use quotes
