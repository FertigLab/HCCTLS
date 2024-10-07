# R version 4.2.2 (2022-10-31 ucrt)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows 10 x64 (build 19042)

#Last revised: 8/9/2024


#.libPaths(paste0(work,"/Libraries"))


rm(list = ls())


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
work<-getwd()

source(paste0(work, "/read_and_cluster_functions.R"))

####REQUIRED LIBRARIES####
library(reshape2)
library(pals)
library(ggplot2)
library(Hmisc)
library(ComplexHeatmap)
library(ggiraphExtra)
library(diffcyt)
library(ggvoronoi)
library(ggpubr)
library(multcomp)
library(sf)
library(clusterSim)
library(circlize)
library(RColorBrewer)
library(stringr)
library(igraph)
library(readxl)
library(dplyr)
library(packcircles)
library(gridExtra)
library(limma)
library(qgraph)
library(flowCore)
library(plot3D)
library(akima)
library(basetheme)
library(pheatmap)

library(ggpubr)
library(ggsci)
library(patchwork)
library(tidyverse)
library(umap)


#======================
####RUNNING DATA####
#======================

####DATA LOADING####

## If loading from prior run, run next line and skip creation of output object at line 73 and 192, and clustering at lines 199 and 206
output<-readRDS('backup_output.rds')

## Read in all files and region selection (for TLS)

tls_md<-read.csv(paste0(work,"/Config/annotated_TLS_with_area_and_type.csv"))
filenames<-paste0(work,"/Data/tlsselection/",tls_md$tls_csv)
csvlist<-lapply(filenames,read.csv)
names(csvlist)<-tls_md$tls_id
tls_ids<-bind_rows(csvlist,.id="tlsid")
tls_ids$sample_id <- substr(tls_ids$tlsid,0,nchar(tls_ids$tlsid)-2)

output <- returnfcs(metaDataFile = paste0(work,"/Config/TLS_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/TLS_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"),
                    tls_md = tls_md,
                    tls_ids = tls_ids)

rm(tls_ids,tls_md) #cleanup

## Set up levels
samplevels=c("TLS_J17136_P02_1", "TLS_J17136_P02_2", 
             "TLS_J17136_P03_1", "TLS_J17136_P03_2", "TLS_J17136_P03_3", "TLS_J17136_P03_4", 
             "TLS_J17136_P07_1", 
             "TLS_J17136_P08A_1", "TLS_J17136_P08B_1", "TLS_J17136_P08B_2", 
             "TLS_J17136_P12_1", "TLS_J17136_P12_2", "TLS_J17136_P12_2", "TLS_J17136_P12_3", 
             "TLS_J17136_P09_1", "TLS_J17136_P09_1", "TLS_J17136_P09_2", "TLS_J17136_P09_3", "TLS_J17136_P09_3", "TLS_J17136_P09_4", "TLS_J17136_P09_4", 
             "TLS_OT1_1", "TLS_OT1_2", "TLS_OT1_3", "TLS_OT1_4", 
             "TLS_OT6A_1", "TLS_OT6A_2", "TLS_OT6A_3", "TLS_OT6A_4", "TLS_OT6A_5", "TLS_OT6A_5", 
             "TLS_OT7_1", "TLS_OT7_2", "TLS_OT7_3", "TLS_OT7_3", "TLS_OT7_4", "TLS_OT7_4", "TLS_OT7_5")

locationlevels = c("Intratumoral", "Peritumoral", "Adjacent Normal")


caselevels=c("02","03","07","08A","08B","12","09","1","6A","7")

tlslevels <- levels(factor(output$tls_md$tls_id))

tls_typelevels=c("mature", "involuted")
typelevels = tls_typelevels

my_comparisons <- list( c("mature", "involuted")) 


####DIAGNOSTICS####
## Spot check - number of cells per sample
cell_table <- table(output$sample_ids)
ggdf <- data.frame(sample_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$case <- factor(output$meta_data$case[match(ggdf$sample_id,output$meta_data$sample_id)], levels = caselevels)
ggdf$type <- factor(output$meta_data$type[match(ggdf$sample_id,output$meta_data$sample_id)], levels= tls_typelevels)
ggdf$location <- factor(output$meta_data$location[match(ggdf$sample_id,output$meta_data$sample_id)], levels=locationlevels)

ggp<-ggplot(ggdf, aes(x = sample_id, y = cell_counts, fill = sample_id)) + 
  geom_bar(stat = 'identity', show.legend = F) + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts.pdf',width=8, height=6);ggp; dev.off()

cell_table <- table(output$tls_ids)
ggdf <- data.frame(tls_id = names(cell_table), 
                   cell_counts = as.numeric(cell_table))
ggdf$case <- factor(output$tls_md$case[match(ggdf$tls_id,output$tls_md$tls_id)], levels = caselevels)
ggdf$type <- factor(output$tls_md$tls_type[match(ggdf$tls_id,output$tls_md$tls_id)], levels= tls_typelevels)
ggdf$location <- factor(output$tls_md$location[match(ggdf$tls_id,output$tls_md$tls_id)], levels=locationlevels)

ggp<-ggplot(ggdf, aes(x = tls_id, y = cell_counts, fill = tls_id)) + 
  geom_bar(stat = 'identity', show.legend = F) + 
  #geom_text(aes(label = cell_counts), angle = 45, hjust = 0.5, vjust = -0.5, size = 2) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size=5)) +  
  scale_x_discrete(drop = FALSE)
pdf('Diagnostics_cellcounts_tls.pdf',width=8, height=6);ggp; dev.off()

## Multi-dimensional scaling plot to show similarities between samples
## Get the mean marker expression per sample
expr_mean_sample_tbl <- data.frame(sample_id = output$tls_ids, fsApply(output$fcstls,exprs)) %>%
  group_by(sample_id) %>%  summarize_all(funs(mean))
expr_mean_sample <- t(expr_mean_sample_tbl[, -1])
colnames(expr_mean_sample) <- expr_mean_sample_tbl$sample_id
mds <- plotMDS(expr_mean_sample, plot = FALSE)
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y,
                   sample_id = colnames(expr_mean_sample))
ggdf$case <- factor(output$tls_md$case[match(ggdf$sample_id,output$tls_md$tls_id)], levels = caselevels)
ggdf$type <- factor(output$tls_md$tls_type[match(ggdf$sample_id,output$tls_md$tls_id)], levels = tls_typelevels)
ggdf$location <- factor(output$tls_md$location[match(ggdf$sample_id,output$tls_md$tls_id)], levels = locationlevels)
ggp<-ggplot(ggdf, aes(x = MDS1, y = MDS2, color = case, shape = type)) +
  geom_point(size = 2.5) +
  # geom_text(aes(label = sample_id), check_overlap = T) +
  theme_bw()+
  theme(plot.background = element_rect(fill="white"),
        panel.background = element_rect(fill="white"),
        panel.grid = element_blank(),
        axis.line = element_line(color="black"),
        axis.text = element_text(color="black"),
        axis.title = element_text(color="black"),
        legend.key = element_rect(fill="white"),
        legend.text = element_text(color="black"))
pdf('Diagnostics_MDS_sample_tls_labeled.pdf',width=6, height=6);ggp; dev.off()

## Colors for the heatmap
color <- rev(brewer.rdbu(100))
pdf('Diagnostics_Heatmap_tls.pdf',width=8, height=8)
pheatmap(expr_mean_sample_tbl[,output$cluster_by], 
         labels_row = expr_mean_sample_tbl$sample_id,
         color = color, 
         display_numbers = F,
         number_color = "black", 
         fontsize_number = 3,
         cellwidth = 8,
         cellheight = 9,
         border_color = NA,
         legend = F,
         treeheight_row = 20,
         clustering_method = "average")
dev.off()


####CLUSTERING####

##Revised loading depending on the diagnostics if needed

tls_md<-read.csv(paste0(work,"/Config/annotated_TLS_with_area_and_type.csv"))
filenames<-paste0(work,"/Data/tlsselection/",tls_md$tls_csv)
csvlist<-lapply(filenames,read.csv)
names(csvlist)<-tls_md$tls_id
tls_ids<-bind_rows(csvlist,.id="tlsid")
tls_ids$sample_id <- substr(tls_ids$tlsid,0,nchar(tls_ids$tlsid)-2)

output <- returnfcs(metaDataFile = paste0(work,"/Config/TLS_metadata.xlsx"),
                    panelDataFile = paste0(work,"/Config/TLS_panel.xlsx"),
                    dataDirectory = paste0(work,"/Data"),
                    tls_md = tls_md,
                    tls_ids = tls_ids)
##Clustering

output[(length(output)+1):(length(output)+3)] <- clusterfcs(fcs=output$fcstls, numclusters=50, #50, #40
                                                            scaleoption = T) 
#output$fcs uses arcsin transformed data for entire ROIs
#output$fcstls uses arcsin transformed data but only for tls regions
#output$fcs2 uses scaled expression data that is scaled for each flowFrame for entire ROIs
#scaleoption scales the dataset across the channels so that the channels with the highest intensities will not dominate the clustering

names(output)[(length(output)-2):(length(output))] <- c('code_clustering','cell_clustering','metaclusters')

####ANNOTATIONS AND VISUALIZATION OF CLUSTERS####

## Load merge file
## Assign colors
clusterMergeFile = paste0(work,"/Config/TLS_merge.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile)

clusterlevels=c(1:50) #40)

colorassigned<-kovesi.rainbow_bgyrm_35_85_c69(length(clusterlevels))

clusternames<-clusterlevels
names(colorassigned)<-clusternames
mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m

## Metacluster heatmaps
plot_clustering_heatmap_wrapper(fcs=output$fcstls,
                                number_clusters = length(unique(output$cell_clustering1m)),
                                cell_clustering = output$cell_clustering, 
                                cluster_by = output$cluster_by,
                                clusterMergeFile = clusterMergeFile,
                                rngmin=0.01, rngmax=0.99,
                                fileName = 'Clusteringheatmap_fcs_tls_50.pdf'); dev.off()


# Check the location of the clusters
exprtb<-as.data.frame(fsApply(output$fcstlsraw,exprs))
exprtb$sample_ids<-tls_ids$sample_id
exprtb$cell_clustering<-output$cell_clustering1m
# tls_ids = tls_ids %>% group_by(sample_id) %>% arrange(.,sample_id,Parameter_1) #the exprtb output has the cell ids in order grouped by fcs rather than TLS, so here we have to rearrange the tls_ids so it's in the same order 
 if (all(unique(exprtb$sample_ids) == unique(tls_ids$sample_id)) == T){ #this doublechecks that we did this correctly, should return all trues
   exprtb$tls_ids<-tls_ids$tlsid
   print("OK to proceed. tls_ids$sample_id is in correct order to match order of samples in exprtb")
   } else {print("tls_ids$sample_id is not in correct order to match order of samples in exprtb")} #note this if_else will not catch if the parameter 1 cell ids are not in the correct order
exprtb$tls_type<-factor(output$tls_md$tls_type[match(exprtb$tls_ids,output$tls_md$tls_id)], levels = tls_typelevels)

##add column to exprtb with original mcd viewer ID
tls_df <- read_excel(path=paste0(work,"/Data/Updated Catalogue of TLS with S numbers.xlsx"))
exprtb$patient_id_ROI = tls_df$Name[match(exprtb$tls_ids, tls_df$ImageJ_ROI)]
rm(tls_df)

## note this object, exprtb, has a column sample_ids which is not in the correct order since it is organized by sample. the column with the correct order is tls_ids, which is arranged by sample/ROI 
pdf("Clusterloccheck_initial.pdf",width=12,height=12)
for(i in 1:length(clusterlevels)){
  clustersubset<-exprtb[exprtb$cell_clustering==i,]
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=tls_type))+
    ggtitle(clusterlevels[i])+
    geom_point(size=0.1)+
    theme_bw()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    facet_wrap(~tls_ids)
  print(ggp)
  };dev.off()

pdf("Clusterloccheck_by_sample_ID_50.pdf",width=12,height=12)
for(i in 1:length(tlslevels)){
  clustersubset<-exprtb[exprtb$tls_ids==tlslevels[i],]
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=tls_type))+
    ggtitle(paste0(tlslevels[i], " | ", unique(clustersubset$patient_id_ROI)," | ", unique(clustersubset$tls_type)))+
    geom_point(size=0.1)+
    theme_bw()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    facet_wrap(~cell_clustering)
  print(ggp)
};dev.off()

## Once annotated, input the revised merge file
clusterMergeFile = paste0(work,"/Config/TLS_merge_DS.xlsx") #create dummy merger numbers prior to annotation
cluster_merging <- read_excel(clusterMergeFile, n_max = 50) #
cluster_merging =  cluster_merging %>% select(original_cluster, new_cluster)

level = 2 #set level of granularity, options are 1 (simplified) or 2 (more granular)
if (level == 1){
    clusterlevels=c("B", "Plasma cell","Th", "Tc","GZMB+ T cell","Proliferating T and B","Mature DC","Macrophage","Stroma","HEV","Tumor")
  } else if (level == 2) {
    clusterlevels = sort(unique(cluster_merging$new_cluster))

    clusterlevels=c(
    "B_BCL6high","B_BCL6low","B_AID+", "Plasma cell",
    "Th_CXCR3low",
    "Tph",
    "Th_CD57+",
    "Th_FOXP3+",
    "GZMB+ T cell","Proliferating T and B",
    "Tc",
    "Macrophage",
    "Mature DC",
    "Stroma","HEV","Tumor", "Unassigned"
    )
    } else print("Please set the identity level before proceeding")

colorassigned<-pals::glasbey(length(unique(cluster_merging$new_cluster)))

clusternames<-clusterlevels
names(colorassigned)<-clusternames
colorassigned["Unassigned"]<-"light gray"
colorassigned %>% pals::pal.bands()

mm1 <- match(output$cell_clustering, cluster_merging$original_cluster)
cell_clustering1m <- cluster_merging$new_cluster[mm1]
output$cell_clustering1m <- cell_clustering1m
#add column for new cluster
exprtb$new_cluster = cluster_merging$new_cluster[match(exprtb$cell_clustering,cluster_merging$original_cluster)]

################## cleanup for TLS_J17136_P03_3_1 ##########################
#here we have to do a manual reassignment of a subset of cells in cluster 21 in TLS_J17136_P03_3_1, which captures multiple cell types
#note we are changing only cell_clustering1m assignment, which we just gave in the couple lines above. but don't change the cell_clustering assignment (i.e. 1:50)

#first, out of 61393 indices in output$tls_ids, identify which ones belong to TLS_J17136_P03_3_1, n=2200
reassignment_indices = which(output$tls_ids=="TLS_J17136_P03_3_1") 
#find indices of all cells assigned to cluster 21
cell_clustering_21 = which(output$cell_clustering==21) 
#to determine which indices in cellclustering1m need to be replaced, subset cell_clustering_21 for those found in reassignment_indices
J17136_P03_3_1_cluster_21 = cell_clustering_21[cell_clustering_21 %in% reassignment_indices]

output$cell_clustering1m[J17136_P03_3_1_cluster_21]#first check their identity prior to reassignment (should be Th_CXCR3low)
exprtb[exprtb$tls_ids=="TLS_J17136_P03_3_1" & exprtb$cell_clustering==21,]$new_cluster #and make the same change in the exprtb object

replacement_cluster_for_J17136_P03_3_1_level1 = "B"
replacement_cluster_for_J17136_P03_3_1_level2 = "B_BCL6low"

if (level == 1){
  output$cell_clustering1m[J17136_P03_3_1_cluster_21] <- replacement_cluster_for_J17136_P03_3_1_level1 #reassign according to the level
  exprtb[exprtb$tls_ids=="TLS_J17136_P03_3_1" & exprtb$cell_clustering==21,]$new_cluster <- replacement_cluster_for_J17136_P03_3_1_level1 #reassign according to the level 
  # clusterlevels = c(replacement_cluster_for_J17136_P03_3_1_level1,clusterlevels)
} else if (level == 2) {
  output$cell_clustering1m[J17136_P03_3_1_cluster_21] <- replacement_cluster_for_J17136_P03_3_1_level2 #reassign according to the level 
  exprtb[exprtb$tls_ids=="TLS_J17136_P03_3_1" & exprtb$cell_clustering==21,]$new_cluster <- replacement_cluster_for_J17136_P03_3_1_level2 #reassign according to the level 
  # clusterlevels = c(replacement_cluster_for_J17136_P03_3_1_level2,clusterlevels)
} else print("Please set the identity level before proceeding")

#then check that this was successful
output$cell_clustering1m[J17136_P03_3_1_cluster_21]
exprtb[exprtb$tls_ids=="TLS_J17136_P03_3_1" & exprtb$cell_clustering==21,]$new_cluster

################## make plots #######################################
plot_clustering_heatmap_wrapper2(fcs=output$fcstls, 
                                 #DS: I changed above line from fcs = output$fcs to fcs = output$fcstls due to error: Error in data.frame(expr, cell_clustering = cell_clustering, check.names = FALSE) arguments imply differing number of rows: 225050, 90344;
                                 # output$cell_clustering has same # of rows as output$fcstls, not output$fcs
                                 colorbar = kovesi.diverging_bwr_40_95_c42(100),
                                 subtype_markers = c(
                                   "VIM",
                                   "CD45",
                                   "CD45RA","CD45RO",
                                   "HLADR","CD86","CD25", 
                                   "CD21", "CD23", "BCL6",
                                   "CD20", "AID","CD138",
                                   "CD3",
                                   "CD4",  
                                   "CXCR3", "CXCR5",   
                                   "ICOS",
                                   "TOX","PD1",
                                   "CD57","FOXP3",
                                   "GZMB",
                                   "KI67",
                                   "CD8","LAG3","CD137", 
                                   "CD69","CD68", "CD11c",
                                   "DCLAMP","CCR7",
                                 "PDPN",
                                 "SMA",
                                  "PNAd", 
                                  "CK","PDL1"
                                   
                                   ),
                                 color_clusters = colorassigned,
                                 cell_clustering = factor(output$cell_clustering1m, levels=clusterlevels), 
                                 fileName = 'Clusteringheatmap_Tc_merged_fcs.pdf');dev.off()


############################ remake plots from above #######################################

pdf("Clusterloccheck_final.pdf",width=12,height=12)
for(i in 1:length(clusterlevels)){
  clustersubset<-exprtb[exprtb$new_cluster==clusterlevels[i],]
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=tls_type))+
    ggtitle(clusterlevels[i])+
    geom_point(size=0.1)+
    theme_bw()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    facet_wrap(~tls_ids)
  print(ggp)
}; dev.off()

#make plot
pdf("Clusterloccheck_by_sample_ID_final.pdf",width=20,height=12)
for(i in 1:length(tlslevels)){
  clustersubset<-exprtb[exprtb$tls_ids==tlslevels[i],]
  title = paste0(tlslevels[i], " | ", unique(clustersubset$patient_id_ROI)," | ", unique(clustersubset$tls_type))
  clustersubset$old_new = paste0(clustersubset$new_cluster, " | ", clustersubset$cell_clustering)
  
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ #this is the big plot
    ggtitle(title)+
    geom_point(size=2)+
    theme_bw()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    guides(fill=guide_legend(ncol=1))+
    scale_color_manual(values=colorassigned) #,limits=clusterlevels)
  
  ggp1.5<- ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ 
    geom_point(size=1,show.legend=F)+
    theme_bw()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    scale_color_manual(values=colorassigned) #,limits=clusterlevels)
  ggp2 <- ggp1.5 + facet_wrap(~new_cluster)

  ggp3 <- ggp+ggp2+plot_layout(guides="collect") &theme(legend.position = 'right')
  print(ggp3)

  # these two lines print an additional page that splits up the plots by cluster, can comment out to get the pdf to compile faster
  ggp4<- ggp1.5 + ggtitle(title)+facet_wrap(~old_new)
  print(ggp4)
}; dev.off()

pdf("Clusterloccheck_by_sample_ID_final_v2.pdf",width=8,height=8)
for(i in 1:length(tlslevels)){
  clustersubset<-exprtb[exprtb$tls_ids==tlslevels[i],]
  title = paste0(tlslevels[i], " | ", unique(clustersubset$patient_id_ROI)," | ", unique(clustersubset$tls_type))
  clustersubset$old_new = paste0(clustersubset$new_cluster, " | ", clustersubset$cell_clustering)
  
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ #this is the big plot
    ggtitle(str_to_sentence(unique(clustersubset$tls_type))
            )+
    geom_point(size=2)+
    # theme_bw()+
    ggprism::theme_prism()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    guides(fill=guide_legend(ncol=1))+
    scale_color_manual(values=colorassigned) #,limits=clusterlevels)
  
  print(ggp)+
    theme(text=element_text(size=16),legend.title = element_blank(),legend.text=element_text(size=16))
  
  # these two lines print an additional page that splits up the plots by cluster, can comment out to get the pdf to compile faster
}; dev.off()

tlslevels_subset = tlslevels[tlslevels%in% c("TLS_OT1_3_1","TLS_OT7_4_1")]

pdf("Clusterloccheck_by_sample_ID_final_subset.pdf",width=18,height=8)
for(i in 1:length(tlslevels_subset)){
  clustersubset<-exprtb[exprtb$tls_ids==tlslevels_subset[i],]
  # title = paste0(tlslevels_subset[i], " | ", unique(clustersubset$patient_id_ROI)," | ", unique(clustersubset$tls_type))
  # clustersubset$old_new = paste0(clustersubset$new_cluster, " | ", clustersubset$cell_clustering)
  
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ #this is the big plot
    ggtitle(str_to_sentence(unique(clustersubset$tls_type)))+
    geom_point(size=2.5)+
    # theme_bw()+
    ggprism::theme_prism()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    guides(fill=guide_legend(ncol=1))+
    scale_color_manual(values=colorassigned,limits=clusterlevels)
  
  ggp1.5<- ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ 
    geom_point(size=1,show.legend=F)+
    # theme_bw()+
    ggprism::theme_prism()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    scale_color_manual(values=colorassigned,limits=clusterlevels)
  ggp2 <- ggp1.5 + facet_wrap(~new_cluster)
  
  ggp3 <- ggp+ggp2+plot_layout(guides="collect") &theme(legend.position = 'right',
                                                        legend.text=element_text(size=16))
  print(ggp3)
  
  # these two lines print an additional page that splits up the plots by cluster, can comment out to get the pdf to compile faster
  # ggp4<- ggp1.5 + ggtitle(title)+facet_wrap(~old_new)
  # print(ggp4)
}; dev.off()

pdf("Clusterloccheck_by_sample_ID_final_v2_subset.pdf",width=8,height=8)
for(i in 1:length(tlslevels_subset)){
  clustersubset<-exprtb[exprtb$tls_ids==tlslevels_subset[i],]
  title = paste0(tlslevels_subset[i], " | ", unique(clustersubset$patient_id_ROI)," | ", unique(clustersubset$tls_type))
  clustersubset$old_new = paste0(clustersubset$new_cluster, " | ", clustersubset$cell_clustering)
  
  ggp<-ggplot(clustersubset, aes(x=X_coord,y=Y_coord,color=new_cluster))+ #this is the big plot
    ggtitle(str_to_sentence(unique(clustersubset$tls_type))
    )+
    geom_point(size=2)+
    # theme_bw()+
    ggprism::theme_prism()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    guides(fill=guide_legend(ncol=1))+
    scale_color_manual(values=colorassigned) #,limits=clusterlevels)
  
  print(ggp)+
    theme(text=element_text(size=16),legend.title = element_blank(),legend.text=element_text(size=16))
  
  # these two lines print an additional page that splits up the plots by cluster, can comment out to get the pdf to compile faster
}; dev.off()

pdf("Clusterloccheck_by_sample_ID_final_v3.pdf",width=12,height=6)

  ggp1<-ggplot(exprtb[exprtb$tls_ids%in%c("TLS_OT1_3_1","TLS_OT7_4_1"),], aes(x=X_coord,y=Y_coord,color=new_cluster))+ #this is the big plot
    # ggtitle("Mature")+
    geom_point(size=1.5)+
    # theme_bw()+
    ggprism::theme_prism()+
    scale_y_reverse()+ #puts plot in mcd viewer orientation
    guides(fill=guide_legend(ncol=1))+
    scale_color_manual(values=colorassigned,limits=clusterlevels)+
    facet_wrap(~tls_type,scales="free",ncol=2)+
   theme(text=element_text(size=12),legend.text=element_text(size=12))
  print(ggp1)
  dev.off()

## this code checks the marker expression for all clusters 
pdf("qc_by_cluster.pdf", width=20, height=20)
clusters_for_qc = sort(unique(exprtb$cell_clustering))#[grep("B_|Tc|Th|Proliferating|Panpos",unique(exprtb$new_cluster))]
for (i in 1:length(clusters_for_qc)) {
    df <- exprtb[exprtb$cell_clustering==clusters_for_qc[i],] %>% 
    select(CellID, CD4, CD8, CD20, sample_ids, tls_ids,cell_clustering, new_cluster)
  if (grepl("Tc_", clusters_for_qc[i])==T) {df$index = df$CD8} else{ #make the index column (i.e. order for x axis) CD8 if Tc cluster
    if (grepl("Th_", clusters_for_qc[i])==T) { 
      df$index = df$CD4} else { #make it CD4 if it's a Th cluster
      df$index = df$CD20} #for all other clusters, use CD20
  }
      
  df = df %>% pivot_longer(cols=c("CD4","CD8","CD20"), names_to="marker",values_to="expression")
  df = arrange(df,index)
  df$Cell_Samp = paste0(df$CellID,"_",df$sample_ids) #this new column seems to be necessary to arrange x axis correctly
  df$Cell_Samp = factor(df$Cell_Samp, levels = unique(df$Cell_Samp))
  ggp <- ggplot(df, 
                aes(x=Cell_Samp, y=expression, group=marker,color=marker))+
    geom_line()+
    theme(axis.text.x=element_blank(),
          axis.ticks.x = element_blank())+
    scale_color_aaas()+
    facet_wrap(~tls_ids,scales="free")+
    ggtitle(paste0(clusters_for_qc[i],"_",unique(df$new_cluster)))
  print(ggp)
  rm(ggp)
};dev.off(); rm(clusters_for_qc)


####UMAP####
## Compute umap
#DS 4/13 9:49am, this function doesn't seem to work because the metadata doesn't contain tls information
#to fix this, I am adding another dataframe to output, called output$meta_data_tls, which is output$meta_data merged with the tls list, so there is a distinct row for each tls
meta_data_tls = data.frame(tls_id = unique(output$sample_ids_tls)) 
meta_data_tls$sample_id = substr(meta_data_tls$tls_id,0,nchar(meta_data_tls$tls_id)-2)
output$meta_data_tls = merge(meta_data_tls,by="sample_id",output$meta_data)
rm(meta_data_tls)

umapRes <- do_umap(fcs=output$fcs2, 
                   subtype_markers = output$subtype_markers,
                   sample_ids = output$sample_ids_tls, #I changed this to output$sample_ids to output$sample_ids_tls
                   cell_clustering = output$cell_clustering, 
                   metadata=output$meta_data_tls, #added _tls here too
                   clusterMergeFile=clusterMergeFile,
                   seed = 1234, 
                   ncells=500,
                   sample_subset=NULL)

mm <- match(as.character(umapRes$sample_id), as.character(output[["meta_data_tls"]]$tls_id))
umapRes$sample_id <- factor(output[["meta_data_tls"]]$sample_id[mm], levels=samplevels) 
umapRes$case <- factor(output[["meta_data_tls"]]$case[mm], levels = caselevels)
umapRes$type <- factor(output[["meta_data_tls"]]$type[mm], levels=typelevels)
umapRes$location <- factor(output[["meta_data_tls"]]$location[mm], levels=locationlevels)

##not needed to run if there aren't any "NA" clusters
umapRes<-umapRes[umapRes$cell_clustering!="NA",]

# dev.off()
pdf('Umaps2.pdf',width=8,height=8)
plotUmap(umapRes = umapRes,
         code_clustering=cell_clustering1m,
         color_clusters = colorassigned,
         subtype_markers = output$subtype_markers)
dev.off()

samplestoexclude=c() 
include_inds <- which(output$tls_ids %nin% samplestoexclude)

## Save output list
saveRDS(output, file="backup_output.rds")
saveRDS(umapRes, file="backup_umap.rds")
# umapRes<-readRDS('backup_umap.rds')

#export metadata
#add pub.id column to output$meta_data_tls (used in plots below)
output$meta_data_tls$pub.id = output$meta_data_tls$sample_id %>% 
  str_replace("TLS_J17136_|TLS_","") %>% 
  str_replace("A_.*|B_.*|_.*", "")

pub.id.levels = output$meta_data_tls$pub.id %>% unique 
output$meta_data_tls$pub.id = factor(output$meta_data_tls$pub.id, levels = pub.id.levels) 

#create table to export # of TLS examined per patient (for paper extended data table)
summary_for_pub <- output$meta_data_tls %>% select(pub.id, type,location)
summary_for_pub
summary_for_pub$type = factor(summary_for_pub$type,levels=c("mature","involuted"),labels=c("Mature","Involuted"))
table(summary_for_pub$pub.id,summary_for_pub$type) 

t1 <- table(summary_for_pub[summary_for_pub$type=="Mature",]$pub.id,summary_for_pub[summary_for_pub$type=="Mature",]$location) 
t1
t1 %>% as.data.frame.matrix() %>% rownames_to_column(var="Patient") %>% 
  mutate(Total = rowSums(.[2:4])) %>% 
  writexl::write_xlsx("summary_of_TLS_for_pub_mature.xlsx")

t2 <- table(summary_for_pub[summary_for_pub$type=="Involuted",]$pub.id,summary_for_pub[summary_for_pub$type=="Involuted",]$location) 
t2
t2 %>% as.data.frame.matrix() %>% rownames_to_column(var="Patient") %>% 
  # mutate(Total = rowSums(.[2:4])) %>% 
  writexl::write_xlsx("summary_of_TLS_for_pub_involuted.xlsx")

t1_plus_t2 <- full_join(rownames_to_column(as.data.frame.matrix(t1),var="Patient"),
          rownames_to_column(as.data.frame.matrix(t2),var="Patient"),
          by = join_by(Patient),
          suffix=c("_Mature", "_Involuted")) %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) sum(.) else "Total"))
t1_plus_t2
t1_plus_t2 %>% writexl::write_xlsx("summary_of_TLS_for_pub_mature_and_involuted.xlsx")

df = summary_for_pub %>% group_by(pub.id,type,location) %>% 
  summarize(n = n()) %>% 
  pivot_wider(names_from = location,values_from=n) %>% 
  rename(Patient = pub.id, `TLS type` = type) %>% relocate(Peritumoral, .after=Intratumoral) %>% 
  mutate_all(~replace(., is.na(.), 0))  
df
writexl::write_xlsx(df,"summary_of_TLS_for_pub.xlsx")

# ## replaces old location column with location2
# output$meta_data_tls$location = output$meta_data_tls$location2
# locationlevels = location2levels

####ABUNDANCE PLOTS####
## Proportion calculations
counts_table <- table(output$cell_clustering1m[include_inds], output$sample_ids_tls[include_inds])
props_table <- t(t(counts_table) / colSums(counts_table)) * 100
counts <- as.data.frame.matrix(counts_table)
props <- as.data.frame.matrix(props_table)

areas <- read.csv(paste0(work,'/Config/annotated_TLS_with_area_and_type.csv'))
areas$area = areas$area/1e6 #convert square micrometer to square millimeter
areas = areas[match(colnames(counts), areas$tls_id),]  #reorders area df to make sure it matches column order in counts df
if (all(colnames(counts) == areas$tls_id)) {
  densities <- t(t(counts)/areas$area)
} else {print("densities computation aborted due to mismatch in order")}

write.csv(counts,'Results_counts.csv')
write.csv(props,'Results_props.csv')
write.csv(densities, 'Results_densities.csv')

## Set up the data frame for proportional plotting
ggdf <- melt(data.frame(cluster = rownames(props),props, check.names = FALSE),
             id.vars = "cluster", value.name = "proportion", 
             variable.name = "tls_id")
ggdf$cluster <- factor(ggdf$cluster, levels=clusterlevels)
ggdf$tls_id <- factor(ggdf$tls_id,levels=tlslevels)
# ggdf$sample_id <- factor(ggdf$sample_id, levels=samplevels)
ggdf$case <- factor(output$meta_data_tls$case[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels = caselevels)
ggdf$type <- factor(output$meta_data_tls$type[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels=typelevels)
ggdf$location <- factor(output$meta_data_tls$location[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels=locationlevels)
# Set up dataframe for densities
 ggdfd <- melt(data.frame(cluster = rownames(densities),densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities",
              variable.name = "tls_id")
 ggdfd$cluster <- factor(ggdfd$cluster, levels=clusterlevels)
 ggdfd$tls_id <- factor(ggdfd$tls_id, levels=tlslevels)
 ggdfd$case <- factor(output$meta_data_tls$case[match(ggdfd$tls_id,output$meta_data_tls$tls_id)], levels = caselevels)
 ggdfd$type <- factor(output$meta_data_tls$type[match(ggdfd$tls_id,output$meta_data_tls$tls_id)], levels=typelevels)
 ggdfd$location <- factor(output$meta_data_tls$location[match(ggdfd$tls_id,output$meta_data_tls$tls_id)], levels=locationlevels)
 
 ##############% CELLS by type, then location ################
ggp2<-ggplot(data=ggdf[ggdf$cluster != "Unassigned",], aes(x=type, y=proportion, color=type))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        # legend.text = element_text(size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
        )+
    stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = my_comparisons#,
#                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                       )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_aaas()+  facet_wrap(~cluster,scales="free")
ggp2

pdf('Abundance_box2_Tc_proportion_merged.pdf',width=12,height=12);ggp2;dev.off();rm(ggp2)

# densities by TLS type
ggp2<-ggplot(data=ggdfd[ggdfd$cluster != "Unassigned",], aes(x=type, y=densities, color=type))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black",lwd=0.5,fatten=1)+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Density (cells/mm²)")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        # legend.text = element_text(size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = my_comparisons#,
                     #symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_aaas()+  facet_wrap(~cluster,#nrow=3,
                                  scales="free")
pdf('Abundance_box2_Tc_densities_merged.pdf',width=12,height=12); ggp2;dev.off();rm(ggp2)

##### by location #######
# subset analysis for mature intratumoral vs mature peritumoral
my_comparisons_location = list(c("Intratumoral", "Peritumoral"))
loc_colors = c("Mature_Intratumoral" = "green4", "Mature_Peritumoral" = "brown3")
# %cells for location
ggdf_loc = ggdf[ggdf$type == "mature" & ggdf$location %in% c("Intratumoral","Peritumoral"),]
ggdf_loc$location = factor(ggdf_loc$location, levels=c("Intratumoral", "Peritumoral"), labels=c("Mature_Intratumoral", "Mature_Peritumoral"))
ggp2<-ggplot(data=ggdf_loc[ggdf_loc$cluster != "Unassigned",], aes(x=location, y=proportion, color=location))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("% of Cells")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        legend.text = element_text(face="bold"),#size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = list(levels(ggdf_loc$location))#,
                     #                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_manual(values=loc_colors)+  
  facet_wrap(~cluster,scales="free")
ggp2
pdf('Abundance_box2_Tc_proportion_merged_location.pdf',width=12,height=12);ggp2;dev.off();rm(ggp2)

# densities by location
ggdfd_loc = ggdfd[ggdfd$type == "mature" & ggdfd$location %in% c("Intratumoral","Peritumoral"),]
ggdfd_loc$location = factor(ggdfd_loc$location, levels=c("Intratumoral", "Peritumoral"), labels=c("Mature_Intratumoral", "Mature_Peritumoral"))

ggp2<-ggplot(data=ggdfd_loc[ggdfd_loc$cluster != "Unassigned",], aes(x=location, y=densities, color=location))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black",lwd=0.5,fatten=1)+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Density (cells/mm²)")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        legend.text = element_text(face="bold"),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = list(levels(ggdfd_loc$location))#,
                     #symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_manual(values=loc_colors)+  
  facet_wrap(~cluster,#nrow=3,
                                  scales="free")+
  ggtitle("Mature TLS cell cluster density by location")
ggp2
pdf('Abundance_box2_Tc_densities_merged_location.pdf',width=12,height=12); ggp2;dev.off();rm(ggp2)
##############################

 # for densities, but subset
coi = c("B_BCL6high", "B_BCL6low", "B_AID+", "Plasma cell",  "Tph", "Th_FOXP3+", "GZMB+ T cell", "Proliferating T and B", "Tc", "Mature DC", "HEV", "Tumor")
non_coi = unique(ggdf$cluster)[!unique(ggdf$cluster)%in% coi]
non_coi = non_coi[non_coi!="Unassigned"]

ggp3<-ggplot(data=ggdfd[ggdfd$cluster %in% coi,], aes(x=type, y=densities, color=type))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black",lwd=0.5,fatten=1)+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Density (cells/mm²)")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=12, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12,face="bold"),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
        axis.title.y = element_text(size=12,face="bold"),axis.text.y=element_text(size=10,face="bold")
  )+
  stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = my_comparisons#,
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_aaas()+  facet_wrap(~cluster,nrow=4,
                                  scales="free")
ggp3
pdf('Abundance_box2_Tc_densities_merged_coi.pdf',width=10,height=12); ggp3;dev.off();rm(ggp3)


ggp4<-ggplot(data=ggdfd[ggdfd$cluster %in% non_coi,], aes(x=type, y=densities, color=type))+
  geom_boxplot(outlier.shape = NA, fill="white",color="black",lwd=0.5,fatten=1)+#,width=0.15,colour = "black" #,fill = "white"
  geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Density (cells/mm²)")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=12, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        legend.text = element_text(size=12,face="bold"),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.adj",vjust=-0.1,comparisons = my_comparisons#,
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_color_aaas()+  facet_wrap(~cluster,nrow=1,
                                  scales="free")
ggp4
pdf('Abundance_box2_Tc_densities_merged_not_coi_for_supplemental.pdf',height=4,width=12); ggp4;dev.off();rm(ggp4)

## Line plots
ggdfl <- melt(data.frame(cluster = rownames(densities),densities, check.names = FALSE),
              id.vars = "cluster", value.name = "densities", 
              variable.name = "tls_id")
ggdfl$cluster <- factor(ggdfl$cluster, levels=clusterlevels)
ggdfl$tls_id <- factor(ggdfl$tls_id, levels=tlslevels)
ggdfl$case <- factor(output$meta_data_tls$case[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels=caselevels)
ggdfl$type <- factor(output$meta_data_tls$type[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels=typelevels)
ggdfl$location <- factor(output$meta_data_tls$location[match(ggdf$tls_id,output$meta_data_tls$tls_id)], levels=locationlevels)

## Set up ggdfl to include liver samples ONLY
# ggdfl_liver<-ggdfl[ggdfl$sample_id %nin% nonliversamples,]

## Line Plot - Mean Case Density
linegraph_data <- ggdfl  %>% group_by(cluster, case, type) %>% summarize_at(vars(densities), funs(mean))
# linegraph_data_pairedonly <-linegraph_data[linegraph_data$case %nin% c("09","18", "21"),]

colorassigned2<-kovesi.rainbow_bgyrm_35_85_c69(length(unique(ggdfl$case)))


lp <- ggplot (data=linegraph_data[linegraph_data$cluster != "Unassigned",], aes(x=type, y = densities, group = case, colour = case)) +
  geom_line(size = 0.75)+
  geom_point(size = 2) +
  facet_wrap(~ cluster, ncol=4,scales = "free") +
  labs(
    x = "TLS type",
    y = "Density (cells/mm²)",
    title = "Cell Clusters, by patient",
    group = "Case"
  ) +
  scale_color_manual(name = "Case", values = colorassigned2)+
  # scale_y_sqrt()
  stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

## Line Plot Aesthetics
lp_styled <- lp + theme_bw() +
  theme(
    plot.title = element_text(face = "bold", size = 14, colour = "black"),
    axis.ticks = element_line(colour = "grey70", size = 0.2),
    panel.grid.major = element_blank(),#element_line(colour = "grey70", size = 0.2),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(), #element_rect(fill = "black", color = "black", size = 1),
    strip.text = element_text(face = "bold", size = 9,color = "black"),
    legend.background = element_rect(fill = "white", colour = "white")
)
lp_styled
pdf('Line_Plot_merged_WH.pdf',width=10,height=10);lp_styled;dev.off();rm(lp);rm(lp_styled)

## Stacked bars
bp <- ggplot(ggdf, aes(x = tls_id, y = proportion, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))

bp.mature <- ggplot(ggdf[ggdf$type=="mature",], aes(x = tls_id, y = proportion, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Mature")+theme(plot.title = element_text(hjust = 0.5))


bp.involuted <- ggplot(ggdf[ggdf$type=="involuted",], aes(x = tls_id, y = proportion, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(),#color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  ylab("")+#"% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))+
  ggtitle("Involuted")+  theme(plot.title = element_text(hjust = 0.5))
bp.2 <- bp.mature +bp.involuted +patchwork::plot_layout(guides="collect")
bp.3 <- bp.mature/bp.involuted +patchwork::plot_layout(guides="collect")

bp
bp.mature
bp.involuted
pdf('Abundance_stackedbar_1.pdf', width=12, height=8); bp; bp.mature;bp.involuted;bp.2;bp.3;dev.off()

bp.mature.2 <- bp.mature + #ggprism::theme_prism()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

bp.involuted.2 <- bp.involuted + #ggprism::theme_prism()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank())



pdf('Abundance_stackedbar_2.pdf', width=12, height=8); 
(bp.mature.2 +bp.involuted.2) +patchwork::plot_layout(guides="collect")

dev.off()

####ABUNDANCE DIFFERENTIAL ANALYSES####

#for comparing response (% cells)

pvalueprops<-data.frame(type=NA,value=NA) #first for AOV

for(i in 1:length(clusterlevels)){
  ggdfsubset<-ggdf[ggdf$cluster==clusterlevels[i],]
  # ggdfsubset<-ggdfsubset[ggdfsubset$response!="UNK",]
  
  res<-aov(proportion ~ type, data = ggdfsubset)
  res<-summary(res)
  res<-unlist(res)["Pr(>F)1"]
  pvalueprops<-rbind(pvalueprops,data.frame(type=clusterlevels[i],value=unlist(res)["Pr(>F)1"]))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")
write.csv(pvalueprops,"Results_props_type_aov.csv")

pvalueprops<-data.frame(type=NA,Var1=NA,Var2=NA,value=NA) #t test for PR vs CR only

for(i in 1:length(clusterlevels)){
  ggdfsubset<-ggdf[ggdf$cluster==clusterlevels[i],]
  # ggdfsubset<-ggdfsubset[ggdfsubset$response %nin% c("SD","UNK"),]
  
  res<-pairwise.t.test(ggdfsubset$proportion, ggdfsubset$type, paired = F, p.adjust.method = "none")
  pvalueprops<-rbind(pvalueprops,cbind(type=rep(clusterlevels[i],nrow(melt(res$p.value))),melt(res$p.value)))
  pvalueprops<-pvalueprops[!is.na(pvalueprops$value),]
}
rownames(pvalueprops)<-1:nrow(pvalueprops)
pvalueprops$FDRadjP<-p.adjust(pvalueprops$value, method="BH")

write.csv(pvalueprops,"Results_props_type_t.csv")

####EXTRACTING MEAN EXPRESSIONS OF KEY MARKERS AND COLORLEGENDS BASED ON EXPRESSION####

#list out the cell types and create a legends data frame

allcelltypes<-clusterlevels
legendctype<-as.data.frame(cbind(paste0("ctype",1:length(allcelltypes)),allcelltypes))
legendctype$maintype<-1
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tc|Th|Tph|T cell")]<-"T cell"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Proliferating T and B")]<-"Other lymphocyte"
legendctype$maintype[str_detect(legendctype$allcelltypes,"B_|Plasma")]<-"B cell"
# legendctype$maintype[str_detect(legendctype$allcelltypes,"NK")]<-"Lym"
legendctype$maintype[str_detect(legendctype$allcelltypes,"DC")]<-"Myeloid"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Macrophage")]<-"Myeloid"
# legendctype$maintype[str_detect(legendctype$allcelltypes,"Gran")]<-"Myl"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Stroma")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"HEV")]<-"Stroma"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Tumor")]<-"Tumor"
legendctype$maintype[str_detect(legendctype$allcelltypes,"Unassigned")]<-"Other"
levels(legendctype$maintype) <- c("B cell","T cell", "Other lymphocyte", "Myeloid", "Stroma","Tumor","Other")   

markerlist = c("CD4","CD8","CD20","FOXP3", "TOX", "PNAd","CXCR3")
# FOR THE NEXT LINE,4/14/2023 DS I changed this to output$fcstls from output$fcs, and sample_id = output$tls.ids rather than  output$sample_id
exprtbl <- data.frame(fsApply(output$fcstls,exprs)[, markerlist], sample_id = output$sample_ids_tls, cluster = output$cell_clustering1m) 
exprtbl<- exprtbl %>% group_by(cluster) %>%   summarise_if(is.numeric, mean, na.rm = TRUE)
  # summarize_all(funs(mean(., na.rm=TRUE)))
# exprtbl<- exprtbl[,1:(1+length(markerlist))]
exprmtx<-as.matrix(exprtbl[,2:(1+length(markerlist))])
rownames(exprmtx)<-unlist(exprtbl[,1])

pdf("Clusterheatmap_focused.pdf", width=8, height=8)
pheatmap(exprmtx, scale="column")
dev.off()

exprdf<-as.data.frame(scale(exprmtx)[legendctype$allcelltypes,]) #create scaled df

exprdf$legendctype<-legendctype$V1 #add legendctype

#create color bars
CD4colors<-num2col(exprdf$CD4, pal=coolwarm(100))
CD8colors<-num2col(exprdf$CD8, pal=coolwarm(100))
CD20colors<-num2col(exprdf$CD20, pal=coolwarm(100))
PNAdcolors<-num2col(exprdf$PNAd, pal=coolwarm(100))
FOXP3colors <- num2col(exprdf$FOXP3, pal=coolwarm(100))
TOXcolors <- num2col(exprdf$TOX, pal=coolwarm(100))

exprdf$CD4colors<-CD4colors
exprdf$CD8colors<-CD8colors
exprdf$CD20colors<-CD20colors
exprdf$PNAdcolors<-PNAdcolors
exprdf$FOXP3colors <- FOXP3colors
exprdf$TOXcolors <- TOXcolors

####DISTANCE RELATIONSHIPS####

#how many of each cell types are there in the dataset?

ggdft <- melt(data.frame(cluster = rownames(counts), counts, check.names = FALSE),
              id.vars = "cluster", value.name = "counts", 
              variable.name = "tls_id")
ggdft$tls_id <- factor(ggdft$tls_id, levels=tlslevels)
ggdft$cluster <- factor(ggdft$cluster, levels=clusterlevels)
ggdft$case <- factor(output$meta_data_tls$case[match(ggdft$tls_id,output$meta_data_tls$tls_id)], levels=caselevels)
ggdft$type <- factor(output$meta_data_tls$type[match(ggdft$tls_id,output$meta_data_tls$tls_id)], levels=typelevels)
ggdft$location <- factor(output$meta_data_tls$location[match(ggdft$tls_id,output$meta_data_tls$tls_id)], levels=locationlevels)

totalcounts<-ggdft %>% group_by(cluster, case, type) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts,"Totalcounts.csv")

totalcounts_tlstype <- ggdft %>% group_by(cluster, type) %>% summarize_at(vars(counts),funs(sum))

write.csv(totalcounts_tlstype, "Totalcounts_tlstype_met.csv")

#percentage of each respective total

totalcounts_mature <- totalcounts[totalcounts$type=="mature",]
totalcounts_involuted <- totalcounts[totalcounts$type=="involuted",] 
totalinvoluted<-sum(totalcounts_involuted$counts)
totalmature<-sum(totalcounts_mature$counts)
pct_involuted <- totalcounts_involuted$counts/totalinvoluted*100
names(pct_involuted)<- totalcounts_involuted$cluster
pct_mature <- totalcounts_mature$counts/totalmature*100
names(pct_mature)<- totalcounts_mature$cluster

pct_mature = tapply(pct_mature, names(pct_mature), sum)
pct_mature = pct_mature[clusterlevels] #reorder following clusterlevels
pct_involuted = tapply(pct_involuted, names(pct_involuted), sum)
pct_involuted = pct_involuted[clusterlevels]

ggdf_pct_involuted <- melt(as.matrix(pct_involuted)) %>% rename(cluster = Var1) %>% select(-Var2)
ggdf_pct_involuted$type<-"involuted"
ggdf_pct_mature <- melt(as.matrix(pct_mature)) %>% rename(cluster = Var1) %>% select(-Var2)
ggdf_pct_mature$type<-"mature"
ggdf_pct<-rbind(ggdf_pct_mature,ggdf_pct_involuted)
# ggdf_pct$cluster<-c(rownames(ggdf_pct_mature),rownames(ggdf_pct_involuted))
rownames(ggdf_pct)<-1:nrow(ggdf_pct)
ggdf_pct$type<-factor(ggdf_pct$type, levels=c("mature","involuted"))

bp <- ggplot(ggdf_pct, aes(x = type, y = value, fill=cluster, order=cluster)) +
  geom_bar(stat = "identity", position="fill", width=0.8) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust = 1, color="black", size=6),
        axis.text.y = element_text(color="black"),
        axis.ticks = element_line(size=0.25),
        strip.text=element_text(size=8),
        strip.background = element_rect(fill=NA, color=NA),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.key.size = unit(.75,'lines'),
        legend.text = element_text(size=8),
        legend.key = element_rect(fill="white"),
  ) +
  xlab("TLS type")+ylab("% of Total Cells")+
  scale_fill_manual(values = colorassigned,
                    breaks = clusterlevels,
                    labels = clusterlevels)+
  scale_y_continuous(expand = c(0,0))+
  guides(fill=guide_legend(ncol=1))
bp
pdf('Abundance_stackedbar_mature_v_involuted.pdf'); bp; dev.off()

# i reused code from comparison of luminal vs triple negative breast cancer, thus need to create LUM and TNC objects
pct_LUM = pct_mature #pulls in pct_mature from above
pct_TNC = pct_involuted

########################################
#get out expression levels, X, and Y coords

expr <- as.data.frame(fsApply(output$fcs1, exprs)) #create expression matrix   #this was won's original code, which contains all 200K+ cells. there is no tls only comparable fcs1 object in output, so i will extract the 90330 cells in tls using join with the columns that are the same in the object exprtb, which we made above
exprtb.subset = exprtb %>% select(CellID,Area,Eccentricity,Perimeter,MajAxis,MinAxis,X_coord,Y_coord, sample_ids,cell_clustering,tls_ids,tls_type)
expr_tls = left_join(exprtb.subset,expr, by=c("CellID","Area","Eccentricity","Perimeter","MajAxis","MinAxis","X_coord","Y_coord"))

expr0<-data.frame(expr_tls[,c(union(output$subtype_markers,output$functional_markers),"CellID","Area","Eccentricity","Perimeter","MajAxis","MinAxis","X_coord","Y_coord","sample_ids","cell_clustering","tls_ids","tls_type"
                          )]#,
                  # cluster=output$cell_clustering1m, #commented this out since I already havea column for clusters, instead will substitute new_cluster names for these numbers
                  # sample_id=output$sample_ids_tls) #commented this out since I already have sample_ids and tls_ids columns in expr_tls
)
expr0$cluster = cluster_merging$new_cluster[match(expr0$cell_clustering,cluster_merging$original_cluster)]
expr0$cluster <- factor(expr0$cluster, levels=clusterlevels)
expr0$tls_ids <- factor(expr0$tls_ids, levels=tlslevels)
expr0$sample_ids <- factor(expr0$sample_ids, levels=samplevels)
expr0$location <- factor(output$meta_data_tls$location[match(expr0$tls_ids,output$meta_data_tls$tls_id)], levels=locationlevels)

expr0<-expr0[expr0$cluster!="NA",]

####################################################
expr0_marker <- pivot_longer(expr0, cols = names(expr0)[1:44], names_to = "marker", values_to = "expression") #make version of expr0 with new columns marker and expression
expr0_marker$marker %>% unique %>% length
expr0_marker = expr0_marker[expr0_marker$marker %in% output$functional_markers,] 
expr0_marker$marker %>% unique %>% length

ggp2<-ggplot(data=expr0_marker, aes(x=tls_type, y=expression))+
  geom_violin(aes(fill=tls_type))+
  geom_boxplot(outlier.shape = NA, width=0.2,fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
  # geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Marker expression")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        # legend.text = element_text(size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_fill_aaas()+
  facet_wrap(~marker,scales="free")+ggtitle("Marker expression across all cells, by TLS type")
# ggp2
pdf('Marker_expr_by_TLS_type_allClusters.pdf',width=12,height=12); ggp2;dev.off();rm(ggp2)

### use lapply to plot marker expression by cluster
ggp2_list <- lapply (levels(expr0_marker$cluster), function(i) {
  ggplot(data=expr0_marker[expr0_marker$cluster==i,], aes(x=tls_type, y=expression))+
  geom_violin(aes(fill=tls_type))+
  geom_boxplot(outlier.shape = NA, width=0.2,fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
  # geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Marker expression")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        # legend.text = element_text(size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
#                     symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                     )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_fill_aaas()+
  facet_wrap(~marker,scales="free")+
    ggtitle(paste0(i, " cluster marker expression"))
})
pdf('Marker_expr_by_TLS_type_byCluster.pdf',width=12,height=12); ggp2_list;dev.off(); rm(ggp2_list)

#plot exhaustion and memory markers per cluster
### for loop to plot marker expression by cluster
expr0_marker_t = expr0_marker[expr0_marker$marker %in% c("CD45RO", "CD25", "CD69", "CD137", "LAG3", "PD1", "TOX"),]
expr0_marker_t$marker = factor(expr0_marker_t$marker, levels=c("CD45RO", "CD25", "CD69", "CD137", "LAG3", "PD1", "TOX"))

ggp2_list <- lapply (levels(expr0_marker_t$cluster), function(i) {
  ggplot(data=expr0_marker_t[expr0_marker_t$cluster==i,], aes(x=tls_type, y=expression))+
    geom_violin(aes(fill=tls_type))+
    geom_boxplot(outlier.shape = NA, width=0.2,fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
    # geom_jitter(width=0.2, size=1.5)+
    xlab("TLS type")+ylab("Marker expression")+
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(), #axis.text = element_text(color="black"),
          # strip.background = element_rect(fill="gray"),
          strip.text = element_text(size=9, face="bold",color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          # legend.key.size = unit(2,'lines'),
          # legend.text = element_text(size=8),
          # legend.key = element_rect(fill="white"),
          axis.ticks.x = element_blank(),
    )+
    stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
#                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                       )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
    scale_fill_aaas()+
    facet_wrap(~marker,scales="free",nrow=1)+
    ggtitle(i)
})
names(ggp2_list) = levels(expr0_marker_t$cluster)
pdf('Marker_expr_by_TLS_type_byCluster_t.pdf',width=12,height=4); ggp2_list;dev.off()

markers_mdc = c("CD11c", "CCR7", "DCLAMP", "HLADR", "CD86")
expr0_marker_mdc = expr0_marker[expr0_marker$marker %in% markers_mdc,]
expr0_marker_mdc$marker = factor(expr0_marker_mdc$marker, levels=markers_mdc)
expr0_marker_mdc = expr0_marker_mdc[expr0_marker_mdc$cluster == "Mature DC",]

ggp3<- ggplot(data=expr0_marker_mdc, aes(x=tls_type, y=expression))+
    geom_violin(aes(fill=tls_type))+
    geom_boxplot(outlier.shape = NA, width=0.2,fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
    # geom_jitter(width=0.2, size=1.5)+
    xlab("TLS type")+ylab("Marker expression")+
    theme(axis.text.x = element_blank(),
          panel.grid = element_blank(), #axis.text = element_text(color="black"),
          # strip.background = element_rect(fill="gray"),
          strip.text = element_text(size=9, face="bold",color="black"),
          panel.background = element_rect(fill="white"),
          legend.title = element_blank(),
          # legend.key.size = unit(2,'lines'),
          # legend.text = element_text(size=8),
          # legend.key = element_rect(fill="white"),
          axis.ticks.x = element_blank(),
    )+
    stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
                       #                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
    )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
    scale_fill_aaas()+
    facet_wrap(~marker,scales="free",nrow=1)+
    ggtitle("Mature DC")
print(ggp3)
pdf('Marker_expr_by_TLS_type_byCluster_mdc.pdf',width=12,height=4); ggp3;dev.off()

pdf('Marker_expr_by_TLS_type_byCluster_mdc_t_pub.pdf',width=11,height=11)
((ggp3+xlab("")+ylab("")) / 
  ggp2_list$Tph / 
  (ggp2_list$Tc +xlab("")+ylab(""))
) + patchwork::plot_layout(guides="collect")
dev.off()

pdf('Marker_expr_by_TLS_type_byCluster_tumor.pdf',width=12,height=4); ggp2_list;dev.off()

markers_tumor = c("HLADR")
expr0_marker_tumor = expr0_marker[expr0_marker$marker %in% markers_tumor,]
expr0_marker_tumor$marker = factor(expr0_marker_tumor$marker, levels=markers_tumor)
expr0_marker_tumor = expr0_marker_tumor[expr0_marker_tumor$cluster == "Tumor",]

ggp4<- ggplot(data=expr0_marker_tumor, aes(x=tls_type, y=expression))+
  geom_violin(aes(fill=tls_type))+
  geom_boxplot(outlier.shape = NA, width=0.2,fill="white",color="black")+#,width=0.15,colour = "black" #,fill = "white"
  # geom_jitter(width=0.2, size=1.5)+
  xlab("TLS type")+ylab("Marker expression")+
  theme(axis.text.x = element_blank(),
        panel.grid = element_blank(), #axis.text = element_text(color="black"),
        # strip.background = element_rect(fill="gray"),
        strip.text = element_text(size=9, face="bold",color="black"),
        panel.background = element_rect(fill="white"),
        legend.title = element_blank(),
        # legend.key.size = unit(2,'lines'),
        legend.text = element_text(face="bold"),#size=8),
        # legend.key = element_rect(fill="white"),
        axis.ticks.x = element_blank(),
  )+
  stat_compare_means(method="wilcox.test",label = "p.value",vjust=-0.1,comparisons = my_comparisons#,
                     #                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
  )+
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
  scale_fill_aaas()+
  facet_wrap(~marker,scales="free",nrow=1)+
  ggtitle("Tumor")
print(ggp4)
pdf('Marker_expr_by_TLS_type_byCluster_tumor.pdf',width=4,height=5); ggp4;dev.off()


####################################################

expr0_LUM<-expr0[expr0$tls_type=="mature",] #here will use mature in place of LUM
expr0_TNC<-expr0[expr0$tls_type=="involuted",] #here will use involuted in place of TNC

###CREATE DISTANCE MATRICES FOR PROGRESSION CRITERIA IN LUM or TNC ONLY

LUM<-unique(expr0[expr0$tls_type=="mature",]$tls_ids)
TNC<-unique(expr0[expr0$tls_type=="involuted",]$tls_ids)

##LUM

expr_LUM<-c()

for(k in 1:length(LUM)){
  expr_k<-expr0[expr0$tls_ids==LUM[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-paste0("ctype",1:length(allcelltypes))
  expr_k<-data.frame(expr_k,d=dummy)

  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for(i in 1:nrow(expr_k)){
    #looking at cell i in core k
    celldistances1<-data.frame(dist=dist1[,i],cluster=expr_k$cluster)
    celldistances1$clustercode<-match(celldistances1$cluster,legendctype$allcelltypes)
    for(j in 1:length(allcelltypes)){
      #among all cells of cell type j, what non-zero (non-self) distance is the closest
      min_d_j<-min(celldistances1[celldistances1$cluster==allcelltypes[j],]$dist[celldistances1[celldistances1$cluster==allcelltypes[j],]$dist!=0])
      expr_k[i,paste0("d.ctype",j)]<-min_d_j
    }
  }
  
  expr_k$ctype_no<-paste0("ctype",match(expr_k$cluster,legendctype$allcelltypes))
  
  expr_k_m<-as.matrix(expr_k[,colnames(expr_k)[str_detect(colnames(expr_k),"d.c")]])
  rownames(expr_k_m)<-expr_k$ctype_no
  colnames(expr_k_m)<-legendctype$V1
  expr_k_m[is.infinite(expr_k_m)]<-NA
  
  expr_LUM<-rbind(expr_LUM,expr_k_m)
  
}

saveRDS(expr_LUM,'backup_dist_LUM_mature.rds')

#load previously saved distance matrix

expr_LUM<- readRDS('backup_dist_LUM_mature.rds')

#create expression + distance data frames

expr_LUMcombined <- cbind(expr0_LUM,expr_LUM)

# DS commented this out 
# expr_LUM_CK<-expr_LUMcombined[expr_LUMcombined$cluster %in% output$functional_markers,]

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_LUMrm0<-expr_LUM[,colnames(expr_LUM) %nin% colnames(expr_LUM)[colSums(expr_LUM, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.01%) cell types 
#I commented this out 5/8/2023 DS since it seems to drop some of the cell types I care about (e.g. B_AID+ is present in mature but absent in involuted, which I want to keep)
#instead just use rm0 object above
expr_LUMrm =expr_LUMrm0
# expr_LUMrm<-expr_LUMrm0[rownames(expr_LUMrm0) %nin% exclude_LUM, colnames(expr_LUMrm0) %nin% exclude_LUM]

##TNC

expr_TNC<-c()

for(k in 1:length(TNC)){
  expr_k<-expr0[expr0$tls_ids==TNC[k],] 
  
  #create placer cols for shortest distance to each cell type
  dummy<-matrix(nrow=nrow(expr_k),ncol=length(allcelltypes)) 
  colnames(dummy)<-paste0("ctype",1:length(allcelltypes))
  expr_k<-data.frame(expr_k,d=dummy)
  
  #get distances from index cell (row) to all cells
  dist1<-as.matrix(dist(cbind(expr_k$X_coord,expr_k$Y_coord)))
  dist1<-as.data.frame(dist1)
  
  for(i in 1:nrow(expr_k)){
    #looking at cell i in core k
    celldistances1<-data.frame(dist=dist1[,i],cluster=expr_k$cluster)
    celldistances1$clustercode<-match(celldistances1$cluster,legendctype$allcelltypes)
    for(j in 1:length(allcelltypes)){
      #among all cells of cell type j, what non-zero (non-self) distance is the closest
      min_d_j<-min(celldistances1[celldistances1$cluster==allcelltypes[j],]$dist[celldistances1[celldistances1$cluster==allcelltypes[j],]$dist!=0])
      expr_k[i,paste0("d.ctype",j)]<-min_d_j
    }
  }
  
  expr_k$ctype_no<-paste0("ctype",match(expr_k$cluster,legendctype$allcelltypes))
  
  expr_k_m<-as.matrix(expr_k[,colnames(expr_k)[str_detect(colnames(expr_k),"d.c")]])
  rownames(expr_k_m)<-expr_k$ctype_no
  colnames(expr_k_m)<-legendctype$V1
  expr_k_m[is.infinite(expr_k_m)]<-NA
  
  expr_TNC<-rbind(expr_TNC,expr_k_m)
  
}

saveRDS(expr_TNC,'backup_dist_TNC_involuted.rds')

#load previously saved distance matrix

expr_TNC<-readRDS('backup_dist_TNC_involuted.rds')

#create expression + distance data frames

expr_TNCcombined <- cbind(expr0_TNC,expr_TNC)

# expr_TNC_CK<-expr_TNCcombined[expr_TNCcombined$cluster %in% output$functional_markers,]

#improving robustness of distance relationships in the dataset by removing cell types without relationships or very rare cell types that would be overrepresented (sampling bias)

#remove cell type columns where there are no cells

expr_TNCrm0<-expr_TNC[,colnames(expr_TNC) %nin% colnames(expr_TNC)[colSums(expr_TNC, na.rm = T)==0]]

#remove rows/columns where there are very rare (<0.5%) cell types
#I commented this out 5/8/2023 DS since it seems to drop some of the cell types I care about (e.g. B_AID+ is present in mature but absent in involuted, which I want to keep)
expr_TNCrm =expr_TNCrm0
# expr_TNCrm<-expr_TNCrm0[rownames(expr_TNCrm0) %nin% exclude_TNC, colnames(expr_TNCrm0) %nin% exclude_TNC]

#################################################################
### 3D plots for proximity based on ECAD and VIM 
akimainterp <- function(exprdf,ctype){
  exprdf_f<-exprdf[!is.na(exprdf[[ctype]]),]
  akima::interp(x=exprdf_f$PNAd,exprdf_f$KI67,1/exprdf_f[[ctype]],duplicate="mean")
}

#for perspective 3D plots
int_ctype1<-akimainterp(exprdf=expr_LUMcombined,ctype="ctype1")
int_ctype11<-akimainterp(exprdf=expr_LUMcombined,ctype="ctype11")
int_ctype12<-akimainterp(exprdf=expr_LUMcombined,ctype="ctype12")
int_ctype13<-akimainterp(exprdf=expr_LUMcombined,ctype="ctype13")
int_ctype14<-akimainterp(exprdf=expr_LUMcombined,ctype="ctype14")


###NETWORK VISUALIZATION 
###by a distance object (of matrix) of shortest distance between cell types

#for LUM
mat_LUM=aggregate(x=expr_LUMrm, by=list(rownames(expr_LUMrm)), FUN=mean, na.rm=T)
groupnames<-mat_LUM$Group.1
mat_LUM<-as.matrix(mat_LUM[,2:ncol((mat_LUM))])
rownames(mat_LUM)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_LUM)<-str_sub(colnames(mat_LUM),6,) #simplify colnames
mat_LUMex<-mat_LUM[colnames(mat_LUM)[c(colnames(mat_LUM) %nin% which(legendctype$allcelltypes=="Unassigned"))], 
                   rownames(mat_LUM)[c(rownames(mat_LUM) %nin% which(legendctype$allcelltypes=="Unassigned"))]
                   ] #excluding Other subtypes  #ds 5/8/2023 
dist_LUM<-dist(mat_LUMex) %>% as.matrix() #i needed to add  as.matrix to make this work DS 5/8/2023

#for TNC
mat_TNC=aggregate(x=expr_TNCrm, by=list(rownames(expr_TNCrm)), FUN=mean, na.rm=T)
groupnames<-mat_TNC$Group.1
mat_TNC<-as.matrix(mat_TNC[,2:ncol((mat_TNC))])
rownames(mat_TNC)<-str_sub(groupnames,6,) # add the rownames
colnames(mat_TNC)<-str_sub(colnames(mat_TNC),6,) #simplify colnames
# mat_TNCex<-mat_TNC[colnames(mat_TNC)[c(colnames(mat_TNC) %nin% "17")], rownames(mat_TNC)[c(rownames(mat_TNC) %nin% "17")]] #excluding Unassigned subtypes #ds 5/8/2023 
mat_TNCex<-mat_TNC[colnames(mat_TNC)[c(colnames(mat_TNC) %nin% which(legendctype$allcelltypes=="Unassigned"))], 
                   rownames(mat_TNC)[c(rownames(mat_TNC) %nin% which(legendctype$allcelltypes=="Unassigned"))]
                   ] #excluding Unassigned subtypes #ds 5/8/2023 commented out
dist_TNC<-dist(mat_TNCex)%>% as.matrix() #i needed to add  as.matrix to make this work DS 5/8/2023

#color of legends revised
legendctypeex<-legendctype[legendctype$maintype!="Other",] #excluding Other subtypes
cols <- brewer.pal(nlevels(as.factor(legendctypeex$maintype)), "Set2")
colorlist <- cols[as.numeric(as.factor(legendctypeex$maintype))]

#generate visualizations of distance relationships using network qgraph

dev.off()

pdf("Distance_Relationships_Mature.pdf",width=8,height=8)
# this indented code was inserted 5/52023 by me, from a slack message from Won
xx<-dist_LUM
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]
colorlistLUM<-colorlist[as.numeric(colnames(mat_LUM))] #set color
names(colorlistLUM)<-as.character(colnames(mat_LUM)) #set color reference names
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_LUM[g_node_ind] 
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Network Visualization for Mature TLS", 
          layout='spring', repulsion=10,
          maximum=max(1/dist_LUM[upper.tri(1/dist_LUM)], #these distance matrices create a 0 when the comparison is between the sample and itself, which then becomes infinity when inverted; taking upper.tri since upper and lower tri are identical
                      1/dist_TNC[upper.tri(1/dist_TNC)],
                      na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels=legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(.5,2),
          edge.color="darkgray",# edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])
legend("bottomleft",
       legend=levels(as.factor(legendctype$maintype))[levels(as.factor(legendctype$maintype)) %nin% "Other"],
       col = cols ,
       bty = "n", pch=20 , pt.cex = 2, cex = 1,
       text.col="black" , horiz = F)

dev.off()

pdf("Distance_Relationships_involuted.pdf",width=8,height=8)
xx<-dist_TNC
yy<-xx[as.character(sort(as.numeric(rownames(xx)))),#make sure all rows and cols are ordered
       as.character(sort(as.numeric(colnames(xx))))]
colorlistTNC<-colorlist[as.numeric(colnames(mat_TNC))] #set color
names(colorlistTNC)<-as.character(colnames(mat_TNC)) #set color reference names
g<-qgraph(1/yy, DoNotPlot=T) #save qgraph object
g_node_ind<-as.numeric(g$graphAttributes$Nodes$names)
nodesizes<-pct_TNC[g_node_ind]
nodesizessc<-log(1.5+nodesizes)*3 #scale nodesizes
g<-qgraph(1/yy,
          title="Network Visualization for Involuted TLS", 
          layout='spring', repulsion=10,
          maximum=max(1/dist_LUM[upper.tri(1/dist_LUM)], #these distance matrices create a 0 when the comparison is between the sample and itself, which then becomes infinity when inverted; taking upper.tri since upper and lower tri are identical
                      1/dist_TNC[upper.tri(1/dist_TNC)],
                      na.rm=T), #max distance has to be chosen properly to scale across all plots
          vsize=nodesizessc,
          labels= legendctype[g$graphAttributes$Nodes$names,]$allcelltypes,
          label.scale=F,
          label.cex=1,
          curve=1.75, arrows=T, diag=F,
          node.label.offset=c(0.5,2),
          edge.color="darkgray",# edge.width=1.75, fade=T,
          minimum=1/250,
          borders=F,
          color=colorlistLUM[g$graphAttributes$Nodes$names])

dev.off()

dist_ttest<-function(group1=expr_Prm,
                     group2=expr_NPrm,
                     ctype1="ctype1",
                     ctype2="ctype2"){
  c1toc2_g1 = as.vector(group1[rownames(group1)==ctype1,ctype2])
  c1toc2_g2 = as.vector(group2[rownames(group2)==ctype1,ctype2])
  res<-t.test(x=c1toc2_g1,
              y=c1toc2_g2,
              paired=F)
  medians<-c(median(c1toc2_g1, na.rm=T),median(c1toc2_g2, na.rm=T))
  names(medians)<-c("P_median","NP_median")
  estimate<-res$estimate
  names(estimate)<-c("P_mean","NP_mean")
  diff_med<-medians["P_median"]-medians["NP_median"]
  diff_mean<-estimate["P_mean"]-estimate["NP_mean"]
  combine<-cbind(t(estimate), t(medians), diff_med=diff_med, diff_mean=diff_mean, pvalue=res$p.value)
  rownames(combine)<-paste(ctype1,ctype2,sep="_")
  print(combine)
}

res0<-data.frame(P_median=NA, NP_median=NA, P_mean=NA, NP_mean=NA, diff_med=NA, diff_mean=NA, pvalue=NA)

###VIOLIN PLOTS OF DISTANCES
expr_P = expr_LUM
expr_NP = expr_TNC

df_P<-melt(expr_P)
df_P$type <- "mature"
df_NP<-melt(expr_NP)
df_NP$type <- "involuted"
df_combined<-rbind(df_P,df_NP)
colnames(df_combined)<-c("index","target","distance","type")
df_combined<-df_combined[!is.na(df_combined$distance),]

df_combined_list_by_index <- lapply (1:length(unique(df_combined$index)), function(i) {
  df_combined[df_combined$index==unique(df_combined$index)[i],]
})

names(df_combined_list_by_index) = unique(df_combined$index)

#define plotting function
distance_violinplots <- function(df) {
    index_as_celltype = legendctype[match(unique(df$index), legendctype$V1),'allcelltypes'] #this object fetches the celltype of index
  df$target = legendctype[match(df$target, legendctype$V1),'allcelltypes'] #this replaces target column with celltype names rather than ctype1, ctype2, etc.

  p<-ggplot(data=df, aes(x=type, y=distance, fill=type))+
    ggtitle(paste0(index_as_celltype, " distance from other celltypes"))+
    geom_violin()+  
    geom_boxplot(outlier.shape = NA,width = 0.15,fill = "white", colour = "black")+
    xlab("")+
    theme_bw()+
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(), 
      panel.background = element_rect(fill="white"),
      axis.text = element_text(color="black"),
      strip.text = element_text(size=9, face="bold",color="black"),        
      legend.title = element_blank(),
    )+
    stat_compare_means(method="wilcox.test",label = "p.adj", vjust=-0.05,comparisons = my_comparisons#,
#                       symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf), symbols = c("****", "***", "**", "*", "ns"))
                       )+
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)))+
    scale_fill_aaas()+facet_wrap(~ target, scales="free")
  return(p)
}

pdf('Distance_violinplots.pdf',width=14,height=14)
lapply(1:length(df_combined_list_by_index), function(i) {
  distance_violinplots(df_combined_list_by_index[[i]])
})
dev.off()

##additional QC 05/05/2023
pdf("qc_dotplots05052023.pdf")
exprtb[exprtb$new_cluster == "Proliferating T and B",] %>% ggplot(.,aes(x=tls_type,y=CK,color=tls_ids)) + geom_jitter()+ggtitle("Proliferating T and B")
exprtb[exprtb$new_cluster == "Proliferating T and B",] %>% ggplot(.,aes(x=CD3,y=CK,color=tls_ids)) + geom_point()+ggtitle("Proliferating T and B")
exprtb[exprtb$new_cluster == "Proliferating T and B",] %>% ggplot(.,aes(x=CD4,y=CD8,color=tls_ids)) + geom_point()+ggtitle("Proliferating T and B")

exprtb[exprtb$new_cluster == "Th_CD57+",] %>% ggplot(.,aes(x=CD4,y=CD8,color=tls_ids)) + geom_point()+ggtitle("Th_CD57+")

exprtb[exprtb$new_cluster == "GZMB+ T cell",] %>% ggplot(.,aes(x=CD4,y=CD8,color=tls_ids)) + geom_point()+ggtitle("GZMB+ T cell")

exprtb[exprtb$new_cluster == "B_AID+",] %>% ggplot(.,aes(x=CD20,y=CD4,color=tls_ids)) + geom_point()+ggtitle("B_AID+")
exprtb[exprtb$new_cluster == "B_AID+",] %>% ggplot(.,aes(x=CD20,y=AID,color=tls_ids)) + geom_point()+ggtitle("B_AID+")
exprtb[exprtb$new_cluster == "B_AID+",] %>% ggplot(.,aes(x=CD20,y=CD4,color=tls_ids)) + geom_point()+ggtitle("B_AID+")

exprtb[exprtb$new_cluster == "Th_FOXP3+",] %>% ggplot(.,aes(x=tls_type,y=FOXP3,color=tls_ids)) + geom_jitter()+ggtitle("FOXP3+ T cell")
exprtb[exprtb$new_cluster == "Th_FOXP3+",] %>% ggplot(.,aes(x=TOX,y=FOXP3,color=tls_ids)) + geom_point()+ggtitle("FOXP3+ T cell")
p1 <-exprtb[exprtb$new_cluster == "Th_FOXP3+",] %>% ggplot(.,aes(x=tls_type,y=CD4,color=tls_ids)) + geom_jitter()+ggtitle("FOXP3+ T cell")
p2 <-exprtb[exprtb$new_cluster == "Th_FOXP3+",] %>% ggplot(.,aes(x=tls_type,y=CD8,color=tls_ids)) + geom_jitter()+ggtitle("FOXP3+ T cell")
p1+p2
dev.off()

#nearest neighbor analysis
source("2_nearestneighborcode_DS.R")
