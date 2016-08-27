#########################################################################################################
# Author: Roula Tsonaka (s.tsonaka@lumc.nl)
#         Leiden University Medical Center (LUMC)
#########################################################################################################
# Title: LongRNAseqVoom.R
# Aim: R code to analyse longitudinal RNAseq data using voom transformation from limma.
# Notes: 
#      This code can be used to analyse longitudinal RNAseq data via the lme4 package 
#      after the voom transformation of R Bioconductor package limma.
#      The main function is longRNAseq(.) with arguments:
#       - formula.lmer = a two-sided linear formula object to be used in lmer(.) describing both the 
#                        fixed-effects and random-effects part of the model, 
#                        with the response on the left of a ~ operator and the terms, 
#                        separated by + operators, on the right. 
#                        Regarding the random-effects terms only random-intercepts are used at the moment.
#       - data = dataframe containing the variables named in formula and the RNAseq data in each column per genomic feature.
#       - Time = a character string with the variable name in data for the followup time.
#       - ID = a character string with the variable name in data for the subject indicator.
#       - gene.nams = a character string with the gene names in data.
#       - dmat.voom = design matrix to be used in limma.
#
# Date: 26AUGUST2016
###########################################################################################################


# Load main function
source("./MainFunction.R")

# Load data
data. <- read.table("./Dataset.txt")


#####################
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
library(limma)
library(lme4)
library(lmerTest)


# Differential expression testing
diff.expr <- longRNAseq(formula.lmer = "group * time + (1|id)", data = data., 
                        Time = "time", ID = "id", 
                        gene.nams = paste("X", 1:100, sep = ""), 
                        dmat.voom = model.matrix(~ group, data = data.))
  
# Results
diff.expr$pvals # pvalues per coefficient and gene
diff.expr$coefs # fixed effects estimates per gene
diff.expr$disp # estimates of random effects variance per gene