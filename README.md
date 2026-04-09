A function for repairing admixture of two cells in a single droplet during single cell sequencing, notably during Cell-Cell-seq. Corrects the diluted distribution for genes which can reasonably be assigned to one cell or the other.

For more details on the application of this method and Cell-Cell-seq in general, please read more here: Systematic mapping of emergent transcriptional states in interacting single-cell dyads by Cell-Cell-seq, https://www.biorxiv.org/content/10.64898/2026.02.05.704136v1

# Dependencies:
ccRepair is based on the Seurat pipeline in R. GGplot is used for some data visualization.
```
library(Seurat)
library(ggplot2)
```
# Getting started
Download the script and source it to your R environment, then provide your RDS and desired parameters to the Create_ccRepair_Data() function
```
source("2026.04.09_ccRepair_function_v1.0_script.R")

repaired.ratio.df<-Create_ccRepair_Data( Object_data= Object_data,
                                  ExpName="ExpName", ##Title for output files
                                  working_directory="/u/project/kp1/jlangerm/Projects/Dicarlo_Nanovial/Analysis/", ###Place for output files
                                  cell_discriminator = markers$Name, ## Vector containing cell labels that will distinguish SINGLET populations
                                  cell_discriminator_pair = c("CD45NV-Tcell","NV-PC3"), ### Singlet Discriminator strings for TWO populations
                                  fold_change_baseline = 1.5 , ### Fold change baseline for considering "gene allegiance". Raw FC. Too low and you will have bleed, too high and not many genes will be adjusted
                                  coexpression_ceiling = 0.1 , ### A percent measure between 0 and 1, lower for less sample bleed, but may need to be adjusted if bleed is naturally high
                                  expression_cutoff = 0.5,  ### Denoising cutoff, ask consistent genes to be above an average normalized count threshold
                                  consistency_function = "median", ### median, mean, or percent - function for the consistency cutoff. consistency is (gene$Mean+0.1) / (gene$Stdev +0.1). mean and median take upper bound over avg mean or median consistency
                                  consistency_percent = NA ,## if using percent give a 0-1 value, will take top x % of consistent genes ordered by consistency measure
                                  dyads_specific_sample = T ,### if your experiment design has a sample where you know your dyads can only be there pick True and below options, F if you're not sure where they are
                                  dyads_specific_discriminator = markers$Name, # or NA if dyads_specific_sample=F
                                  dyads_specific_name = "NV-PC3-Tcell", ## or NA if dyads_specific_sample=F
                                  dyad_identity_min1 = 0.05,
                                  dyad_identity_min2 = 0.2,
                                  consistency_gene_percent_threshold = 0.1 ,## Percent of consistent genes that must be expressed in a cell to reliably call it a dyad and adjust its expression, 0-1
                                  write_output_data=T)
```


# How It Works

In order to adjust differential transcripts which were depleted by the sequencing saturation of mixed cells in a dyad, we constructed a repaired cell by gene matrix using ccRepair. Singlet samples were used to construct a baseline set of consistent identity genes, based on pairwise differential gene expression between PC3s and T-Cells. Consistent genes with allegiance to each sample were determined by at least 1.5 fold change and 0.05 pvalue, at least 0.5 normalized counts on average, less than 10% of cells expressing in the comparison sample, and finally were selected by taking the top two quartiles of a consistency score represented by (gene mean)/(gene standard deviation). This filtration process found roughly 100 genes per singlet sample which represented well expressed and differential identity genes for constructing a baseline distribution of cell type expression. An average of these top consistent genes is constructed as a reference score for each singlet population, and then for each dyad cell, although for individual cells drop out is controlled by only considering genes greater than zero in the consistency score. The dyad cell score is compared to both singlet reference score, and the ratio of that comparison is used to adjust all genes greater than 1.5 fold change from the original singlet comparison for both sides of the comparison. Dyads are also filtered for a sufficient level of dual identity gene expression in a direct ratio of identity genes - with a cutoff of 5% of T-cell gene expression and 20% of PC3 gene expression. The adjusted differential gene matrix is remerged with unadjusted genes, and then reintegrated with the remaining non-dyad non-adjusted cells to create a corrected gene by cell matrix. 
