A function for repairing admixture of two cells in a single droplet during single cell sequencing, notably during Cell-Cell-seq. Corrects the diluted distribution for genes which can reasonably be assigned to one cell or the other.

For more details on the application of this method and Cell-Cell-seq in general, please read more here:   
### Systematic mapping of emergent transcriptional states in interacting single-cell dyads by Cell-Cell-seq  
https://www.biorxiv.org/content/10.64898/2026.02.05.704136v1

# Dependencies:
For ease of use, ccRepair uses the Seurat pipeline in R. GGplot is used for some data visualization.
```
library(Seurat)
library(ggplot2)
```
# Getting started
Download the script and source it to your R environment, then provide your RDS and desired parameters to the Create_ccRepair_Data() function
```
source("2026.04.09_ccRepair_function_v1.0_script.R")

your_rds <- readRDS("your_seurat_object.rds")
sample_name_vector <- c("Tcell","CancerCell","Tcell","Mixed,"CancerCell"...)


repaired.ratio.df <- Create_ccRepair_Data( Object_data= your_rds,
                                  ExpName="ExpName", ##Title for output files
                                  working_directory="/u/project/Analysis/", ###Place for output files
                                  cell_discriminator = sample_name_vector, ## Vector containing cell labels that will distinguish SINGLET populations
                                  cell_discriminator_pair = c("Tcell","CancerCell"), ### Singlet Discriminator strings for TWO populations
                                  fold_change_baseline = 1.5 , ### Fold change baseline for considering "gene allegiance". Raw FC. Too low and you will have cell cell bleed, too high and not many genes will be adjusted
                                  coexpression_ceiling = 0.1 , ### A percent measure between 0 and 1, lower for less sample bleed, but may need to be adjusted if bleed is naturally high between similar cells
                                  expression_cutoff = 0.5,  ### Denoising cutoff, ask consistent genes to be above an average normalized count threshold
                                  consistency_function = "median", ### median, mean, or percent - function for the consistency cutoff. consistency is (gene$Mean+0.1) / (gene$Stdev +0.1). mean and median take upper bound over avg mean or median consistency
                                  consistency_percent = NA ,## if using percent for function, provide a 0-1 value, will take top x % of consistent genes ordered by consistency measure
                                  dyads_specific_sample = T ,### if your experiment design has a sample where you know your dyads can only be there pick True and below options, F if you're not sure where they are
                                  dyads_specific_discriminator = sample_name_vector, # or NA if dyads_specific_sample=F
                                  dyads_specific_name = "Mixed", ## or NA if dyads_specific_sample=F
                                  dyad_identity_min1 = 0.05, ## Minimum allowable identity gene expression ratio for cell type 1, used to confirm a dyad prior to adjustment. You may need to run this once to determine the best level empirically
                                  dyad_identity_min2 = 0.2, ## Minimum allowable identity gene expression ratio for cell type 2, used to confirm a dyad prior to adjustment. 
                                  consistency_gene_percent_threshold = 0.1 ,## Percent of consistent genes that must be expressed in a cell to reliably call it a dyad and adjust its expression, a value from 0-1
                                  write_output_data=T) ## This will write the ccRepair cell x gene matrix dataframe, the cell metadata with ccRepair consistency gene data, and a list of the top consistent allegiance genes used
```
The most important parameters are the discriminator parameters, which tell the function which cells to take for singlets and which cells to take for dyad, assuming that is know. Sequencing a cell cell seq experiment alone as a single sample is one way to have a priori knowledge of dyads. If it is not known, ccRepair will attempt to guess which dyads are present based on mixing of identity genes from the whole cell population. 
The default parameters represent a good specific baseline, but for your data they may require empirical determination; check the GeneMeans scatter plot to confirm cell-type specific expression of the labelled consistency genes, and that they are enough of them per side of the comparison.



# Follow up
Take the ccRepair dataframe and rerun analysis on it, the Dyad gene expression will be much stronger, allowing differentiation expression, gene correlation, and even projection tools to separate the Cell-cell effects much easier.
```

test_obj <- CreateSeuratObject(counts = repaired.ratio.df, project = "ccRepaired_data")
test_obj@assays$RNA$data <- test_obj@assays$RNA$counts
test_obj = FindVariableFeatures(test_obj, verbose = F)
all.genes <- rownames(test_obj)
test_obj <- ScaleData(test_obj,  vars.to.regress = c("nFeature_RNA"))
test_obj <- RunPCA(test_obj, features = VariableFeatures(test_obj))
test_obj <- FindNeighbors(test_obj, dims = 1:20)
test_obj <- FindClusters(test_obj, resolution = 0.5)
test_obj = RunUMAP(test_obj, dims = 1:20, verbose = T)


cell.design <- test_obj@meta.data
umap.meta <- as.data.frame(test_obj@reductions$umap@cell.embeddings)
names(umap.meta)<- c("X","Y")
cell.design <- cbind(cell.design, umap.meta)

plot1<-ggplot(cell.design)+
  geom_point(aes(X,Y, color=sample_name_vector))+
  theme_classic()
```
![Results of ccRepair adjustment on UMAP projection](https://github.com/Teneth/ccRepair/blob/main/ccrepair_result.png)

# How It Works

Singlet samples are used to construct a baseline set of consistent identity genes, based on pairwise differential gene expression between PC3s and T-Cells. If you cannot be certain of singlet populations in your data, you must either use your best guess or introduce control cell population data from another available dataset. ccRepair requires baseline singlet populations to determine gene allegiance.

With standard parameters, consistent genes with allegiance to each sample are determined by at least 1.5 fold change and 0.05 pvalue, at least 0.5 normalized counts on average, less than 10% of cells expressing in the comparison sample, and finally are selected by taking the top two quartiles of a consistency score represented by the function f(Gene_c)= (gene mean)/(gene standard deviation). Initial DEG discovery uses Seurat's FindAllMarkers() negative binomial function. 
This filtration process should find 25~200 genes representing well expressed and differential identity genes for constructing a baseline distribution of cell type specific expression. An average of these top consistent genes is constructed as a reference score for each singlet population, and then for each dyad cell on a per cell basis. For individual cells, drop out is controlled by only considering genes with normalized counts greater than zero in the consistency score. 
The dyad cell scores are compared to both singlet reference scores, and the ratio of that comparison is used to adjust all genes greater than 1.5 fold change from the original singlet comparison for both sides of the comparison. This adjustment is based on the expression of consistent genes, NOT on the ratio of cell type genes, as that admixture rate is unknown and possibly individual doublet dependent. ccRepair adjusts allegiance genes only using the ratio of cell type consistent genes to the that cell type's singlet standard.
Dyads are further filtered for a sufficient level of dual identity gene expression in a direct ratio of identity genes. In our paper, we observed T-cells and PC3 transcripted mixed roughly 1:4, so we used a cutoff of 5% of T-cell gene expression and 20% of PC3 gene expression. This is to confirm that only true Dyads with sufficient 2-cell signal are fixed. The adjusted differential gene matrix is remerged with unadjusted genes, and then reintegrated with the remaining non-dyad non-adjusted cells to create a corrected gene by cell matrix. 
