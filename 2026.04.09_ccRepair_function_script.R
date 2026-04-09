


library(Seurat)
library(ggplot2)



Create_ccRepair_Data <- function( Object_data= your_rds,
                                     ExpName="ExpName", ##Title for output files
                                     working_directory="/projectdir/", ###Place for output files
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
                                     write_output_matrix=T){ 
  
  
  
  
  if(length(cell_discriminator_pair)!=2){
    print("incorrect # of discriminator classes, must be 2")
    break
  }
  
  if(length(cell_discriminator)!=ncol(Object_data)){
    print("error - cell_discriminator vector length must match cell number")
    break
  }
  
  ##Setup
  
  dir.create(paste0(working_directory,"/",Sys.Date(),"_", ExpName))
  outdir <- c(paste0(working_directory,"/",Sys.Date(),"_", ExpName,"/"))
  
  
  
  Object_data@meta.data$Name <- cell_discriminator
  
  if(dyads_specific_sample==T){
    Object_data@meta.data$Name2 <- dyads_specific_discriminator
  }
  
  cell.metadata <- as.data.frame(Object_data@meta.data)
  
  
  ######
  ##Diff test genes
  
  
  diff.exp.data <-FindAllMarkers(Object_data[,Object_data@meta.data$Name %in% cell_discriminator_pair], group.by="Name", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  
  diff.exp.data<-diff.exp.data[diff.exp.data$avg_log2FC > log(fold_change_baseline,2)   ,]
  
  
  
  ##Compare to gene tester
  # gene.tester <- read.table("/u/project/kp1/jlangerm/Projects/Dicarlo_Nanovial/Analysis/2026-03-09_NV18_GEO_Submission/NV16.StimulatedTcell_vs_PC3.DiffGenes.txt")
  # gene.tester<-gene.tester[abs(log(gene.tester$Fold_Change,10))> log(1.5,10)   ,]
  # diff.exp.data[diff.exp.data$gene %in% head(gene.tester$Gene,25),]
  ##All pointing the right direction
  
  
  
  
  ## Pick best consistent diff genes
  
  ### Side 1
  gene.tester <- diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[1],]
  top.genes <- gene.tester$gene
  
  cells.of.interest<- rownames(cell.metadata[ cell.metadata$Name %in%  cell_discriminator_pair[1],])
  gene.design <- cbind.data.frame(Gene= rownames(Object_data)[rownames(Object_data) %in% top.genes],
                                  Mean= rowMeans(Object_data[rownames(Object_data) %in% top.genes,colnames(Object_data) %in% cells.of.interest]))

  gene.design$Stdev <- apply(Object_data@assays$RNA$data[rownames(Object_data) %in% top.genes,colnames(Object_data) %in% cells.of.interest] ,1, sd)
  gene.design$Consistency <- (gene.design$Mean+0.1) / (gene.design$Stdev +0.1)

  ##Filter
  gene.design <- gene.design[gene.design$Gene %in% gene.tester[gene.tester$pct.2< coexpression_ceiling,]$gene,]
  gene.design <- gene.design[gene.design$Mean> expression_cutoff,]
     #### median, mean, or percent functions
  if(consistency_function=="median"){
  gene.design <-  gene.design[gene.design$Consistency>   median(gene.design$Consistency),]
  }
  if(consistency_function=="mean"){
    gene.design <-  gene.design[gene.design$Consistency>   mean(gene.design$Consistency),]
  }
  if(consistency_function=="percent"){
    gene.design <- gene.design[with(gene.design, order(-Consistency)),]
    gene.design <-  head(gene.design,  round(consistency_percent*nrow(gene.design),1))
  }
  ###Top genes side 1 selected
  top.consistent.side1.genes <- gene.design$Gene
  
  
  
  
  ### Side 2
  gene.tester <- diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[2],]
  top.genes <- gene.tester$gene
  
  cells.of.interest<- rownames(cell.metadata[ cell.metadata$Name %in%  cell_discriminator_pair[2],])
  gene.design <- cbind.data.frame(Gene= rownames(Object_data)[rownames(Object_data) %in% top.genes],
                                  Mean= rowMeans(Object_data[rownames(Object_data) %in% top.genes,colnames(Object_data) %in% cells.of.interest]))
  
  gene.design$Stdev <- apply(Object_data@assays$RNA$data[rownames(Object_data) %in% top.genes,colnames(Object_data) %in% cells.of.interest] ,1, sd)
  gene.design$Consistency <- (gene.design$Mean+0.1) / (gene.design$Stdev +0.1)
  
  ##Filter
  gene.design <- gene.design[gene.design$Gene %in% gene.tester[gene.tester$pct.2< coexpression_ceiling,]$gene,]
  gene.design <- gene.design[gene.design$Mean> expression_cutoff,]
  #### median, mean, or percent functions
  if(consistency_function=="median"){
    gene.design <-  gene.design[gene.design$Consistency>   median(gene.design$Consistency),]
  }
  if(consistency_function=="mean"){
    gene.design <-  gene.design[gene.design$Consistency>   mean(gene.design$Consistency),]
  }
  if(consistency_function=="percent"){
    gene.design <- gene.design[with(gene.design, order(-Consistency)),]
    gene.design <-  head(gene.design,  round(consistency_percent*nrow(gene.design),1))
  }
  ###Top genes side 2 selected
  top.consistent.side2.genes <- gene.design$Gene
  
  ###########
  ###report
  print(paste0("Using ",length(top.consistent.side1.genes)," genes for side 1 and ",length(top.consistent.side2.genes)," genes for side 2"))
  ##Should you be worried?
  if( length(top.consistent.side1.genes)<25|length(top.consistent.side2.genes)<25 ){
    print("Consistent gene numbers are quite low for one of the sides, consider loosening filters")
  }
  
  if(( length(top.consistent.side1.genes)>500|length(top.consistent.side2.genes)>500)){
    print("Consistent gene numbers are quite high for one of the sides, consider tightening filters")
  }



 #####
  ### Create adjusted cell matrix for all cells
  side1_set <-diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[1],]$gene
  side2_set <-diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[2],]$gene


  
  
  ############ Non-Zero Mean model- only use genes for ratio that are detected to adjust for single cell drop out vs all genes detected average comparison
  plot.data <- cell.metadata
  simulcount.df <- Object_data@assays$RNA$data[row.names(Object_data) %in% c(top.consistent.side1.genes), ]
  simulcount.df[simulcount.df>0]<- 1
  plot.data$Score_Side1_Consistent <-  colSums(Object_data@assays$RNA$data[row.names(Object_data) %in% top.consistent.side1.genes,])/colSums(simulcount.df)
  simulcount.df <- Object_data@assays$RNA$data[row.names(Object_data) %in% c(top.consistent.side2.genes), ]
  simulcount.df[simulcount.df>0]<- 1
  plot.data$Score_Side2_Consistent <-  colSums(Object_data@assays$RNA$data[row.names(Object_data) %in% top.consistent.side2.genes,])/colSums(simulcount.df)
  baseline.side1<- mean(plot.data[plot.data$Name==cell_discriminator_pair[1],]$Score_Side1_Consistent)
  baseline.side2<- mean(plot.data[plot.data$Name==cell_discriminator_pair[2],]$Score_Side2_Consistent)
  plot.data$Consistent_Ratio_Side1 <- plot.data$Score_Side1_Consistent/baseline.side1
  plot.data$Consistent_Ratio_Side2 <- plot.data$Score_Side2_Consistent/baseline.side2
  ######################
  
  
  matrix.sided.side1 <- Object_data@assays$RNA$data[rownames(Object_data) %in% side1_set, ]
  matrix.sided.side2 <- Object_data@assays$RNA$data[rownames(Object_data) %in% side2_set, ]
  
  matrix.sided.side1<- t(t(matrix.sided.side1)*(1/plot.data$Consistent_Ratio_Side1))
  matrix.sided.side2<- t(t(matrix.sided.side2)*(1/plot.data$Consistent_Ratio_Side2))
  
  matrix.sided.join <- rbind(matrix.sided.side1 ,matrix.sided.side2 )
  
  
  
  ### Display genes chosen as sample mean vs sample mean
  mean.cell.data <- cbind.data.frame(Gene= rownames(Object_data)[rownames(Object_data) %in% diff.exp.data$gene],
                                     Sample1_Mean= rowMeans(Object_data[rownames(Object_data) %in% diff.exp.data$gene,colnames(Object_data) %in% rownames(cell.metadata[ cell.metadata$Name %in%  cell_discriminator_pair[1],])]),
                                     Sample2_Mean= rowMeans(Object_data[rownames(Object_data) %in% diff.exp.data$gene,colnames(Object_data) %in% rownames(cell.metadata[ cell.metadata$Name %in%  cell_discriminator_pair[2],])]))
  
  plot1 <- ggplot(mean.cell.data)+
    geom_point(aes(Sample1_Mean, Sample2_Mean))+
    geom_point(data=mean.cell.data[mean.cell.data$Gene %in% top.consistent.side1.genes,] ,aes(Sample1_Mean, Sample2_Mean), color="#0a83cf")+
    geom_point(data=mean.cell.data[mean.cell.data$Gene %in% top.consistent.side2.genes,] ,aes(Sample1_Mean, Sample2_Mean), color="green4")+
    # geom_point(data=gene.tester[gene.tester$Gene %in% top.consistent.tcell.genes2,] ,aes(Sample1_Mean, Sample2_Mean), color="blue")+
    labs(x="Side1 Expression", y="Side2 Expression")+
    theme_classic()
  ggsave(plot=plot1, paste0(outdir, ExpName, ".Side1_vs_Side2.ThresholdCon.GeneMeans.scatter.png"),
         height=4, width=5)
  
  
  
  
  if(is.null(Object_data@reductions$pca)==T){
    if( is.null(Object_data@reductions$umap)==T){
      # print("PCA and UMAP found")
      # } else {
      print("UMAP and PCA not found, running UMAP derivation")
      Object_data = FindVariableFeatures(Object_data, verbose = F)
      all.genes <- rownames(Object_data)
      Object_data <- ScaleData(Object_data,  vars.to.regress = c("nFeature_RNA"))
      
      Object_data <- RunPCA(Object_data, features = VariableFeatures(object = Object_data))
      Object_data <- FindNeighbors(Object_data, dims = 1:20)
      Object_data <- FindClusters(Object_data, resolution = 0.5)
      Object_data = RunUMAP(Object_data, dims = 1:20, verbose = T)
      
    }
  }
  
  umap.meta <- as.data.frame(Object_data@reductions$umap@cell.embeddings)
  
  names(umap.meta)<- c("X","Y")
  plot.data$X <- umap.meta$X
  plot.data$Y <- umap.meta$Y
  
  
  ###Find likely Dyads
  ### 
  ##Have to use the top private genes to determine cell score
  test.df <- plot.data
  test.df[is.na(test.df)]<-0
  test.df$Score_side1_private <- colMeans(Object_data@assays$RNA$data[rownames(Object_data) %in% head(diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[1],]$gene,25),])
  test.df$Score_side2_private <- colMeans(Object_data@assays$RNA$data[rownames(Object_data) %in% head(diff.exp.data[diff.exp.data$cluster==cell_discriminator_pair[2],]$gene,25),])
  
  test.df$Score1<- test.df$Score_side1_private/(test.df$Score_side1_private + test.df$Score_side2_private)
  test.df$Score2<- test.df$Score_side2_private/(test.df$Score_side1_private + test.df$Score_side2_private)

  
  # head(test.df[with(test.df, order(Score)),])
  
  # 
  # plot1 <- ggplot(test.df)+
  #   geom_point(data=markers,aes(X,Y), color="grey85")+
  #   geom_point(aes(X,Y, color=Score1))+
  #   scale_color_gradientn(colors=c("grey55", "yellow","orange","red"))+
  #   theme_classic()
  # ggsave(plot=plot1, paste0(outdir, Exp.Name, ".score1_in_NV_PC3.UMAP.png"),
  #        height=4, width=5)
  
  # cells.wanted<- test.df
  # cells.wanted$Score1<- cells.wanted$Score_PC3_Consistent/(cells.wanted$Score_PC3_Consistent + cells.wanted$Score_TCell_Consistent)
  # cells.wanted$Score2<- cells.wanted$Score_TCell_Consistent/(cells.wanted$Score_PC3_Consistent + cells.wanted$Score_TCell_Consistent)
  
  
  
  # dyads_specific_discriminator <- markers$Name # or NA
  # dyads_specific_name <- "NV-PC3-Tcell"
  
  
  
  if(dyads_specific_sample==T){
    cells.wanted<- test.df[test.df$Name2==dyads_specific_name,]
  } else {
    cells.wanted<- test.df
  }
  
  
  # mean(cells.wanted$Score1)
  
  print(paste0("Identity ratio in putative dyads is Score1/Score2= ", round(mean(cells.wanted$Score1)/mean(cells.wanted$Score2),2) ))
 
  if(mean(cells.wanted$Score1)/mean(cells.wanted$Score2)>1.5|mean(cells.wanted$Score1)/mean(cells.wanted$Score2)<0.66){
    print(" Ratio far from 1, suggesting unequal mixing. You should adjust the dyad_identity_min scores to reflect this skew. Does not affect calculation but may add too many non-dyads to adjusted matrix")
  }

  # 
  # 
  # if(mean(cells.wanted$Score1)/mean(cells.wanted$Score2)>1.5){
  #   dyad_identity_min2 <-dyad_identity_min1*mean(cells.wanted$Score1)/mean(cells.wanted$Score2)
  # }
  # 
  # if(mean(cells.wanted$Score1)/mean(cells.wanted$Score2)<0.66)){
  #   dyad_identity_min2 <-dyad_identity_min1*mean(cells.wanted$Score1)/mean(cells.wanted$Score2)
  # }
  # 
  
  # cells.wanted2 <- cells.wanted[cells.wanted$Name=="NV-PC3-Tcell",]
  cells.wanted <- cells.wanted[cells.wanted$Score1> dyad_identity_min1 & cells.wanted$Score2> dyad_identity_min2, ]
  # cells.wanted2 <- cells.wanted2[cells.wanted2$Name=="NV-PC3-Tcell",]
  
  if(nrow(cells.wanted)<10){
    print("Error - Upstream filters have removed too many cells. Stopping")
    break
  }

  
  
  plot1 <- ggplot(cells.wanted)+
    geom_point(data=plot.data,aes(X,Y), color="grey85")+
    geom_point(aes(X,Y), color="red4")+
    # scale_color_gradientn(colors=c("grey55", "yellow","orange","red"))+
    theme_classic()
  ggsave(plot=plot1, paste0(outdir, Exp.Name, ".putative.dyad.mixed.cells.UMAP.png"),
         height=4, width=5)
  
  
  
  
  
  
  
  
  ### Filter cells a little to make sure they are not sneaking in
  
  ###Check that at least X% of consistency genes are expressed
  simulcount.df <- Object_data@assays$RNA$data[row.names(Object_data) %in% c(top.consistent.side1.genes,top.consistent.side2.genes), colnames(Object_data) %in% rownames(cells.wanted)]
  simulcount.df[simulcount.df>0]<- 1

  cells.wanted$Side1_Simulcount <- colSums(simulcount.df[row.names(simulcount.df) %in% c(top.consistent.side1.genes),])
  cells.wanted$Side2_Simulcount <- colSums(simulcount.df[row.names(simulcount.df) %in% c(top.consistent.side2.genes),])
  
 
  
  ##Original filtration criteria but I reasoned using the altered ratio was flawed, need to use the original ratio
  ##   which was already done for the file read into here
  cells.wanted2 <- cells.wanted[cells.wanted$Side1_Simulcount>(consistency_gene_percent_threshold * length(top.consistent.side1.genes)) & cells.wanted$Side2_Simulcount>( consistency_gene_percent_threshold * length(top.consistent.side2.genes)),]

  
  print(paste0("Final filters have identified ", nrow(cells.wanted2), " mixed dyad cells"))
  
  
  matrix.sided.mixed <- matrix.sided.join[,colnames(matrix.sided.join) %in% rownames(cells.wanted2) ]
  
  ## Rejoin to main data
  genes.not.affected <- rownames(Object_data)[!rownames(Object_data) %in% row.names(matrix.sided.mixed)]
  extra.gene.df <-Object_data@assays$RNA$data[rownames(Object_data) %in% genes.not.affected, colnames(Object_data) %in% colnames(matrix.sided.mixed)]
  
  ratioed.df <- rbind(extra.gene.df, matrix.sided.mixed)
  ratioed.df <- ratioed.df[rownames(Object_data),]
  
  cells.not.affected <-colnames(Object_data)[!colnames(Object_data) %in% colnames(matrix.sided.mixed)]
  extra.cell.df <-Object_data@assays$RNA$data[, colnames(Object_data) %in% cells.not.affected]
  
  repaired.ratio.df <- cbind(extra.cell.df,ratioed.df)
  repaired.ratio.df<- repaired.ratio.df[,colnames(Object_data)]
  
  
  if(write_output_matrix==T){
    write.table(repaired.ratio.df, gzfile(paste0(outdir, ExpName, ".ccRepaired.full.dataframe.txt.gz")),
                sep="\t", quote=F)
  }
  
  repaired.ratio.df
}
  
  




























