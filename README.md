# LongRNAseq

Modelling longitudinal RNAseq data in lme4 after voom transformation

* `Dataset.txt`: 
Contains a toy dataset with longitudinal RNAseq collected on 20 subjects at 4 occasions: baseline and at 0.2, 0.4 and 0.6 years, afterwards. 
Subjects are assumed to be assigned to 2 treatments groups A and B and our interest is to test group differences at each time point but also over time.
The data have been simulated from an overdispersed poisson mixed effects model and the genes: 
29 79 41 86 91 are assumed to be differentially expressed between the 2 groups at baseline and over time.

* `MainFunction.R`: Main function to model longitudinal RNAsed data via lme4 after using the voom transformation of limma. 
		  The voom transformation is applied on RNAseq data per time point and then the derived weights are plugged in lmer to test for differential gene expression.
		  The output contains a list with elements: (1) gene names, (2) p-values per coefficient and gene, (3) fixed effects estimates per gene
		  and (4) estimates of random effects variance per gene. In the current version of this fnction only random intercepts are supported to model the within subject correlations.

* `LongRNAseqVoom.R`: Example to model longitudinal RNAsed data via lme4 after using the voom transformation of limma.

