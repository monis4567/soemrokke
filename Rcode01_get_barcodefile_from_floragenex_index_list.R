#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-

#____________________________________________________________________________#
# R-code provided for the project:
# GBS Population genetics on thorny skates
library(readxl)
# define working directory
wd00 <- "/home/hal9000/Documents/shrfldubuntu18/soemrokke"
# define input file
inf01 <- "20220113_WWF-Denmark_Skate_20210727-02025.xlsx"
# paste together path and input flie
pthinf01 <- paste0(wd00,"/",inf01)
# read in excel file as tibble, get specific sheet and skip 20 rows
tibl_inx01 <- readxl::read_xlsx(pthinf01, sheet="Index List", skip=20)
# drop first row
tibl_inx01 <- tibl_inx01[-1,]
#substitute spaces underscores
colnames(tibl_inx01) <- gsub(" ","_",colnames(tibl_inx01) )

# with the aim getting an index barcode list as described under section 
# '4.1.2 Specifying the barcodes'
#https://catchenlab.life.illinois.edu/stacks/manual/#prun
outfl1 = "part03B_barcode_index_list.txt"
# paste together path and input flie
pthoutf01 <- paste0(wd00,"/",outfl1)

# use tab as separator
write.table(tibl_inx01, file=pthoutf01, sep="\t",
            row.names = F, # do not use row names
            col.names = F, # do not use columns names
            quote = F) # do not use quotes

#tail(tibl_inx01)




#