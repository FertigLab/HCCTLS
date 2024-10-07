
####TOP 1-3 NEIGHBORS FROM CK+ CLUSTERS####

#extract expression and neighbors
# here was won's original code, commented out

# exprall <- fsApply(output$fcs1, exprs)
# exprall <- data.frame(exprall[,c("CK8","CK14","VIM","ECAD","CellId","Num_Neighbors","NN1","NN2","NN3")],
#                       objtype=output$cell_clustering1m,
#                       sample_id=output$sample_ids)

# since i need to use the tls specific cells ~61k, i reuse an object
# i reacted ~980, expr_tls, in my version of Won's IMC pipeline script
library(ComplexHeatmap)
exprall = expr_tls 
exprall$cluster = cluster_merging$new_cluster[match(exprall$cell_clustering,cluster_merging$original_cluster)]
exprall$cluster <- factor(exprall$cluster, levels=clusterlevels)
exprall$tls_ids <- factor(exprall$tls_ids, levels=tlslevels)
exprall$sample_ids <- factor(exprall$sample_ids, levels=samplevels)
exprall$location <- factor(output$meta_data_tls$location[match(exprall$tls_ids,output$meta_data_tls$tls_id)], levels=locationlevels)
exprall$tls_type<-factor(output$tls_md$tls_type[match(exprtb$tls_ids,output$tls_md$tls_id)], levels = tls_typelevels)

output$functional_markers
exprall <- data.frame(exprall[,c(
                                # "CD4","CD8","CD20", "CXCR3", "TOX","DCLAMP","CCR7","CK", "FOXP3","GZMB",
                                output$functional_markers,"CellID","N_N","NN_1","NN_2","NN_3","cluster","tls_ids","sample_ids","tls_type")],
                      objtype=output$cell_clustering1m #, 
                      # sample_id=output$sample_ids # I commented out since i already have this in dataframe
)


exprall = rename(exprall, CellId = CellID,Num_Neighbors=N_N,NN1=NN_1,NN2=NN_2,NN3=NN_3) #renames the columns from my dataframe to the names used in Won's code

#note from this line onward, all instanes of exprall$tls_id were originally exprall$sample_id in Won's code. I have batch replaced them for my TLS analysis
rownames(exprall)<-paste(exprall$tls_id,exprall$CellId,sep="_")

exprall$NN1 <- paste(exprall$tls_id,exprall$NN1,sep="_")
exprall$NN2 <- paste(exprall$tls_id,exprall$NN2,sep="_")
exprall$NN3 <- paste(exprall$tls_id,exprall$NN3,sep="_")

NN1ind <- match(exprall$NN1,rownames(exprall))
exprall$N1type <- exprall[NN1ind,]$objtype

NN2ind <- match(exprall$NN2,rownames(exprall))
exprall$N2type <- exprall[NN2ind,]$objtype

NN3ind <- match(exprall$NN3,rownames(exprall))
exprall$N3type <- exprall[NN3ind,]$objtype

exprall<-exprall[exprall$objtype!="NA",]

#choose the cell types desired for visualizing the neighboring relationships to, from CK+ cells
#double check to change the pdf names below
# immunecells<-c("Tc","Th_FOXP3+","Th_CD57+","Th_CXCR3low","GZMB+ T cell","Mature DC","B_I","B_II", "B_III")
##subset for mature type first

#subsetting 
NN_mature <- exprall[exprall$tls_type=="mature",]

#make the following replacements
# CK8 -> HEV
# CK8_E -> Tc
# CK8_V -> Mature DC
# CK8_14_E -> GZMB+ T cell
# CK8_14_EV -> Proliferating T and B

#create vector of clusters to evaluate
clusters.of.interest = unique(exprall$objtype) 

#remove "B and T" cluster
clusters.of.interest = clusters.of.interest[clusters.of.interest%nin%c("Tumor", "Stroma", "Unassigned")]

  #c("HEV", "Tc", "Mature DC", "GZMB+ T cell", "Proliferating T and B")
coi.names = str_replace_all(clusters.of.interest,c(" "="_","\\+"= "","/"="_vs_")) #replaces spaces and + signs 

#define a nearest neighbor dataframe creating function, this is implemented for 3 nearest neighbors
create.nndf <- function (df, clusters.of.interest, coi.names) {
  #first nearest neighbor dataframe
  #create list of first nearest neighbor for every 10k cells (normalized by total first)
  
  mature.list = lapply(1:length(clusters.of.interest), function(i) {
    table(df[df$objtype==clusters.of.interest[i],]$N1type)[clusterlevels] / nrow(df[df$objtype==clusters.of.interest[i],]) * 10000
  }) %>% `names<-`(.,coi.names)
  df1 = lapply(mature.list, as.vector) %>% as.data.frame() #compile data into df1 object
  rownames(df1) <- clusterlevels
  # df1<-df1[immunecells,]
  # rownames(df1) <- paste("N1",rownames(df1),sep="_")
  rownames(df1) <- paste(rownames(df1),"N1", sep="_")
  df1[is.na(df1)]<-0
  
  #second nearest neighbor dataframe
  mature.list = lapply(1:length(clusters.of.interest), function(i) {
    table(df[df$objtype==clusters.of.interest[i],]$N2type)[clusterlevels] / nrow(df[df$objtype==clusters.of.interest[i],]) * 10000
  }) %>% `names<-`(.,coi.names)
  df2 = lapply(mature.list, as.vector) %>% as.data.frame() #compile data into df2 object
  rownames(df2) <- clusterlevels
  # df2<-df2[immunecells,]
  # rownames(df2) <- paste("N2",rownames(df2),sep="_")
  rownames(df2) <- paste(rownames(df2),"N2",sep="_")
  df2[is.na(df2)]<-0
  
  # #third nearest neighbor dataframe
  # mature.list = lapply(1:length(clusters.of.interest), function(i) {
  #   table(df[df$objtype==clusters.of.interest[i],]$N3type)[clusterlevels] / nrow(df[df$objtype==clusters.of.interest[i],]) * 10000
  # }) %>% `names<-`(.,coi.names)
  # df3 = lapply(mature.list, as.vector) %>% as.data.frame() #compile data into df2 object
  # rownames(df3) <- clusterlevels
  # # df3<-df3[immunecells,]
  # rownames(df3) <- paste(rownames(df3),"N3",sep="_")
  # df3[is.na(df3)]<-0
  
  #collate dataframes for all top 3 neighbors
  df<-rbind(df1,df2#,df3
            )
  return(df)
}

NNdf <- create.nndf(NN_mature, clusters.of.interest = clusters.of.interest, coi.names = coi.names)
#clean up low count data (anything with less than 1000 count in the subset can be removed)- primary has plenty
# 
# mfilterselect<-function(objtype=objtype, 
#                         minnum=1000, 
#                         tissuetypes=tissuetypes){
#   tofilter<-names(which(table(NN_maturem[NN_maturem$objtype==objtype,]$site)<minnum))[names(which(table(NN_maturem[NN_maturem$objtype==objtype,]$site)<1000)) %in% tissuetypes]
#   ifelse(tofilter!=0,return(paste(objtype,"met",tofilter,sep="_")),return("nothing_to_filter"))
#   #return(paste(objtype,"met",tofilter,sep="_"))
# }

# totest<-list("HEV","Mature DC","Mature DC","GZMB+ T cell","Proliferating T and B","Proliferating T and B")
# lowcountsexclude<-c()
# for(i in 1:length(totest)){
#   lowcountsexclude<-c(lowcountsexclude,mfilterselect(totest[i],1000,c("LUNG","LIVER","BRAIN")))}
# 
# colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]
# 
# NNdf<-NNdf[,colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]]

#annotations


#dont think I need the annotation column for specific markers since I'm not looking at relative expression of specific markers
markerstoannotate <- c() #c("CD4","CD8","DCLAMP","GZMB")
anno_expr<-NN_mature %>% group_by(objtype) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr<-anno_expr[anno_expr$objtype %in% clusters.of.interest,]
anno_expr$type<- "mature"
anno_expr$name<- anno_expr$objtype %>%  #paste(anno_expr$objtype,"mature",sep="_")
  str_replace_all(c(" "="_","\\+"= "","/"="_vs_")) #this is used to make the format of the rownames in the colannos object below the same as the NNdf column names above
anno_expr = anno_expr[match(names(NNdf), anno_expr$name),] #rearranges the rows according to order of columns in NN_df (see https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-vector-with-specific-order)

# DCLAMP<-anno_expr$DCLAMP
# GZMB<-anno_expr$GZMB
# CD4<-anno_expr$CD4
# CD8<-anno_expr$CD8
type<-anno_expr$type
names<-anno_expr$name

colannos<-data.frame(
                    # CD4=CD4,
                    #  CD8=CD8,
                    #  DCLAMP=DCLAMP,
                    #  GZMB=GZMB,
                     type=type,
                     row.names=sort(names,decreasing=F))
colannos_mature <- colannos

colannos_colors_type <- c("#3B4992E5","#EE0000E5")
names(colannos_colors_type) <- c("mature","involuted")
# colannos_colors_e <- num2col(colannos$CD8, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_e) <- colannos$CD8
# colannos_colors_v <- num2col(colannos$CD4, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_v) <- colannos$CD4
# colannos_colors_8 <- num2col(colannos$DCLAMP, colorRampPalette(c("black","green"))(100), ref=c(0,6.5)) #note range set
# names(colannos_colors_8) <- colannos$DCLAMP
# colannos_colors_14 <- num2col(colannos$GZMB, colorRampPalette(c("black","cyan"))(100), ref=c(0,5.2)) #note range set
# names(colannos_colors_14) <- colannos$GZMB
colannos_colors <- list(type=colannos_colors_type#,
                        # CD8=sort(colannos_colors_e), 
                        # CD4=sort(colannos_colors_v),
                        # DCLAMP=sort(colannos_colors_8),
                        # GZMB=sort(colannos_colors_14)
                        )

#draw heatmap
NNdf_mature<-NNdf[rowSums(NNdf)!=0,] #get rid of rows with all zeros
new_order = sort(colnames(NNdf_mature),decreasing=F) #reorder columns alphabetically
NNdf_mature = NNdf_mature[, new_order]
NNdf_mature = NNdf_mature[sort(rownames(NNdf_mature),decreasing=F),] #sort rownames

#this is commented out since my markerstoannotate is null; if I am including specific markers I need this line
# colannos_mature<-colannos_mature[rownames(colannos_mature) %in% colnames(NNdf_mature),]

#remove unassigned neighbors
row_names_to_remove <- row.names(NNdf_mature)[grep("Unassigned|Epithelial|Stroma",row.names(NNdf_mature))]
NNdf_mature = NNdf_mature[!(row.names(NNdf_mature) %in% row_names_to_remove),]

pdf("NN_mature.pdf")
# df = NNdf_mature
p0<-ComplexHeatmap::pheatmap(t(NNdf_mature), 
         main="Mature",
         name = " ",
         cluster_rows = F, cluster_cols = F, scale = "row", 
         # annotation_row = colannos_mature, 
         # annotation_colors = colannos_colors,
         # cutree_rows = 4,
         color = rev(brewer.rdbu(100)),
         border_color = NA,
         # treeheight_row = 15,
         # treeheight_col = 15,
         # cellheight = 8, 
         cellwidth = 8,
         show_colnames = T)
print(p0)
dev.off()

####repeat for involuted####

#subset involuted
NN_involuted <- exprall[exprall$tls_type=="involuted",]

NNdf <- create.nndf(NN_involuted, clusters.of.interest = clusters.of.interest, coi.names = coi.names)
# #clean up low count data (anything with less than 1000 count in the subset can be removed)- primary has plenty
# 
# mfilterselect<-function(objtype=objtype, 
#                         minnum=1000, 
#                         tissuetypes=tissuetypes){
#   tofilter<-names(which(table(NN_maturem[NN_maturem$objtype==objtype,]$type)<minnum))[names(which(table(NN_maturem[NN_maturem$objtype==objtype,]$type)<1000)) %in% tissuetypes]
#   ifelse(tofilter!=0,return(paste(objtype,"met",tofilter,sep="_")),return("nothing_to_filter"))
#   #return(paste(objtype,"met",tofilter,sep="_"))
# }
# 
# totest<-list("CK14","CK14_E","CK14_EV","DCLAMP")
# lowcountsexclude<-c()
# for(i in 1:length(totest)){
#   lowcountsexclude<-c(lowcountsexclude,mfilterselect(totest[i],1000,c("LUNG","BRAIN")))}
# 
# colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]
# 
# NNdf<-NNdf[,colnames(NNdf)[colnames(NNdf) %nin% lowcountsexclude]]
# 

#annotations

# markerstoannotate <- c("CD4","CD8","DCLAMP","CK14") #instead of using this line, we just use previously defined object above

anno_expr<-NN_involuted %>% group_by(objtype) %>% summarise_at(vars(markerstoannotate), mean)
anno_expr<-anno_expr[anno_expr$objtype %in% clusters.of.interest,]
anno_expr$type<- "involuted"
anno_expr$name<- anno_expr$objtype %>%  #paste(anno_expr$objtype,"mature",sep="_")
  str_replace_all(c(" "="_","\\+"= "","/"="_vs_")) #this is used to make the format of the rownames in the colannos object below the same as the NNdf column names above
anno_expr = anno_expr[match(names(NNdf), anno_expr$name),] #rearranges the rows according to order of columns in NN_df (see https://stackoverflow.com/questions/11977102/order-data-frame-rows-according-to-vector-with-specific-order)

# DCLAMP<-anno_expr$DCLAMP
# GZMB<-anno_expr$GZMB
# CD4<-anno_expr$CD4
# CD8<-anno_expr$CD8
type<-anno_expr$type
names<-anno_expr$name

colannos<-data.frame(
                    # CD4=CD4,
                    #  CD8=CD8,
                    #  DCLAMP=DCLAMP,
                    #  GZMB=GZMB,
                     type=type,
                     row.names=sort(names,decreasing=F))
colannos_involuted<-colannos

# this is commented out because I will use the same colannos_colors object as above
# colannos_colors_p <- c("black","lightgray")
# names(colannos_colors_p) <- c("PRI","MET")
# colannos_colors_s <- brewer.set1(5)[c(1:3,5)]
# names(colannos_colors_s) <- c("mature","MET","LUNG","BRAIN")
# colannos_colors_e <- num2col(colannos$CD8, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_e) <- colannos$CD8
# colannos_colors_v <- num2col(colannos$CD4, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_v) <- colannos$CD4
# colannos_colors_8 <- num2col(colannos$DCLAMP, colorRampPalette(c("black","green"))(100), ref=c(0,6.5)) #note range set
# names(colannos_colors_8) <- colannos$DCLAMP
# colannos_colors_14 <- num2col(colannos$ck14, colorRampPalette(c("black","cyan"))(100), ref=c(0,5.2)) #note range set
# names(colannos_colors_14) <- colannos$ck14
# colannos_colors <- list(met=colannos_colors_p, 
#                         type=colannos_colors_s,
#                         CD8=sort(colannos_colors_e), 
#                         CD4=sort(colannos_colors_v),
#                         DCLAMP=sort(colannos_colors_8),
#                         ck14=sort(colannos_colors_14))

#draw heatmap
NNdf_involuted<-NNdf[rowSums(NNdf)!=0,] #get rid of rows with all zeros
new_order = sort(colnames(NNdf_involuted),decreasing=F) #reorder columns alphabetically
NNdf_involuted = NNdf_involuted[, new_order]
NNdf_involuted = NNdf_involuted[sort(rownames(NNdf_involuted),decreasing=F),] #sort rownames

#this is commented out since my markerstoannotate is null; if I am including specific markers I need this line
# colannos_involuted<-colannos_involuted[rownames(colannos_involuted) %in% colnames(NNdf_involuted),]

#remove unassigned neighbors

row_names_to_remove <- row.names(NNdf_involuted)[grep("Unassigned|Epithelial|Stroma",row.names(NNdf_involuted))]
NNdf_involuted = NNdf_involuted[!(row.names(NNdf_involuted) %in% row_names_to_remove),]

pdf("NN_involuted.pdf")
p1<-ComplexHeatmap::pheatmap(t(NNdf_involuted), 
         main="Involuted",
         name = " ",
         cluster_rows = F, cluster_cols = F, scale = "row", 
         # annotation_row = colannos_involuted, 
         # annotation_colors = colannos_colors,
         # cutree_rows = 4,
         color = rev(brewer.rdbu(100)),
         border_color = NA,
         # treeheight_row = 15,
         # treeheight_col = 15,
         # cellheight = 8, 
         cellwidth = 8,
         show_colnames = T)
print(p1)
dev.off()

p3<-p0+p1
pdf("NN_mature_and_involuted.pdf")
print(p3)
dev.off()
###########################
###JOINT HEATMAP

colnames(NNdf_mature)<-paste0(colnames(NNdf_mature),"_mature")
colnames(NNdf_involuted)<-paste0(colnames(NNdf_involuted),"_involuted")
NNdf_all <- cbind(NNdf_mature,NNdf_involuted)
new_order = sort(colnames(NNdf_all),decreasing=F)
NNdf_all <- NNdf_all[, new_order]

colannos_mature$type<-factor("mature")
rownames(colannos_mature)<-paste0(rownames(colannos_mature),"_mature")
colannos_involuted$type<-factor("involuted")
rownames(colannos_involuted)<-paste0(rownames(colannos_involuted),"_involuted")
colannos_all<-rbind(colannos_mature,colannos_involuted)

#this line removed since markerstoannotate is null
# colannos_all<-colannos_all[,c("CD4","CD8","DCLAMP","GZMB","type")] #rearrange the order in which the annos show

# commented out since this duplicates previous code
# colannos_colors_p <- c("black","lightgray")
# names(colannos_colors_p) <- c("PRI","MET")
# colannos_colors_s <- brewer.set1(5)
# names(colannos_colors_s) <- c("mature","MET","LUNG","LIVER","BRAIN")
# colannos_colors_e <- num2col(colannos_all$CD8, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_e) <- colannos_all$CD8
# colannos_colors_v <- num2col(colannos_all$CD4, colorRampPalette(c("black","yellow"))(100))
# names(colannos_colors_v) <- colannos_all$CD4
# colannos_colors_8 <- num2col(colannos_all$DCLAMP, colorRampPalette(c("black","green"))(100)) 
# names(colannos_colors_8) <- colannos_all$DCLAMP
# colannos_colors_14 <- num2col(colannos_all$ck14, colorRampPalette(c("black","cyan"))(100)) 
# names(colannos_colors_14) <- colannos_all$ck14
# colannos_colors_type <- c("black","lightgray")
# names(colannos_colors_type)<-c("mature","involuted")

# colannos_colors_joint <- list(met=colannos_colors_p, 
#                               type=colannos_colors_s,
#                               type=colannos_colors_type,
#                               CD8=sort(colannos_colors_e), 
#                               CD4=sort(colannos_colors_v),
#                               DCLAMP=sort(colannos_colors_8),
#                               ck14=sort(colannos_colors_14))

NNdf_all = NNdf_all[sort(rownames(NNdf_all),decreasing=F),]

df = t(NNdf_all)
rownames(df)
rbind_result  <- c(rbind( #the rownames are sorted alphabetically, which puts involuted above mature, so this reorders them
  seq(2,nrow(df), by=2),
  seq(1,nrow(df)-1,by=2)))
df = df[rbind_result,]

pdf("NN_all.pdf")
p3 <- pheatmap(df, 
         main="", #Top neighbors of mature and involuted TLS",
         cluster_rows = F, cluster_cols = F, scale = "row", 
         annotation_row = colannos_all,
         annotation_colors = colannos_colors,
         color = rev(brewer.rdbu(100)),
         # cutree_rows = 4,
         border_color = NA,
         # treeheight_row = 5,
         # treeheight_col = 5,
         # cellheight = 8, 
         # cellwidth = 10,
         show_colnames = T,
         gaps_row = c(seq(2, ncol(NNdf_all)-2,by=2))
         # gaps_col = c(seq(2, nrow(NNdf_all)-2,by=2))
         )
print(p3)
dev.off()
