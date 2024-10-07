######################################
######################################
#functions for single patients

tc <- function(x) {
  #set x axis parameters
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  
  #make plots
  tempnames <- names(x$data)
  plot1<-lapply(1:length(tempnames), function(j) {
    trackClonotypes(x$data,list(tempnames[j], 10), .col="aa") %>% vis(., .order = order(x$meta$Order))+
      my_xaxis+
      ggsci::scale_fill_futurama()+xlab("")+
      ggtitle(paste0("Top 10 ", if_else(tcrbcr == "TCRB", paste0("TCR","\U03B2"), "IGH"),
                     " clonotypes for patient ", i)
      )
    })
  
  pdf(paste0(output.path,"tc.pdf"), width=10, height=7)
  print(plot1)
  dev.off()
}

######################
exploratory.single.patient <- function(x) {
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  p1 <- repExplore(x$data, .method = "clones") #number of T cells
  plot1 <- vis(p1)+my_xaxis+  theme(legend.position="NULL")
  
  exp_len <- repExplore(x$data, .method = "len", .col = "aa")
  exp_len_1 <- vis(exp_len) +theme_prism()
  
  exp_cnt <- repExplore(x$data, .method = "count")
  exp_cnt_1 <- vis(exp_cnt) +theme_prism()
  
  exp_vol <- repExplore(x$data, .method = "volume")
  exp_vol_1 <- vis(exp_vol)+my_xaxis+theme(legend.position="NULL")
  temp <- list(plot1, exp_vol_1, exp_len_1,exp_cnt_1)
  
  pdf(paste0(output.path, "exploratory_single_patient.pdf"))
  print(temp)
  dev.off()
}

# clonal.prop method computes the proportion of repertoire occupied by the pools of cell clones
clonality.single.patient <- function(x) {
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  
  imm_pr <- repClonality(x$data, .method = "clonal.prop")
  imm_pr_1 <- vis(imm_pr)+my_xaxis+theme(legend.position="NULL")
  #The top method considers the most abundant cell clonotypes:
  imm_top <- repClonality(x$data, .method = "top", .head = c(10, 100, 1000, 3000, 10000))
  imm_top_1 <- vis(imm_top)+my_xaxis#+theme(legend.position="NULL")
  
  # While the rare method deals with the least prolific clonotypes:
  imm_rare <- repClonality(x$data, .method = "rare")
  imm_rare_1 <- vis(imm_rare)+my_xaxis#+theme(legend.position="NULL")
  
  #Finally, the homeo method assesses the clonal space homeostasis, i.e., the proportion of the repertoire occupied by the clones of a given size:
  imm_hom <- repClonality(x$data, .method = "homeo", .clone.types = c(Small = .0001, Medium = .001, Large = .01, Hyperexpanded = 1))
  imm_hom_1 <- vis(imm_hom)+
    theme_prism()+my_xaxis#+theme(legend.position="NULL")
  
  temp <- list(imm_pr_1,imm_top_1,imm_rare_1,imm_hom_1)
  
  pdf(paste0(output.path, "clonality.pdf"))
  print(temp)
  dev.off()
}


### Repertoire overlap and public clonotypes
repOv.single.patient <- function(x) {
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  
  #create sc, which provides a standard fill scale between 0 and 1 (for morisita plots)
  myPalette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
  sc <- scale_fill_gradientn(colours = myPalette(100), limits=c(0, 1))
  
  
  imm_ov1 <- repOverlap(x$data, .method = "public", .verbose = F)
  plot1 <- vis(imm_ov1, .text.size=4)+
    my_xaxis+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample))
  
  #Public method repertoire overlap by patient
  imm_ov2 <- repOverlap(x$data, .method = "morisita", .col="aa", .verbose = F)
  imm_ov2_1 <- vis(imm_ov2, .text.size = 4)
  plot2 <- vis(imm_ov2, .text.size = 4)+
    my_xaxis+
    sc+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample))
  
  imm_ov2.nt <- repOverlap(x$data, .method = "morisita", .col="nt")
  plot2.nt <- vis(imm_ov2.nt, .text.size = 2)+
    my_xaxis+
    sc+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample)) +ggtitle("Morisita overlap by CDR3nt")
  
  imm_ov2.v <- repOverlap(x$data, .method = "morisita", .col="v")
  plot2.v<- vis(imm_ov2.v, .text.size = 2)+
    my_xaxis+
    sc+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample)) +ggtitle("Morisita overlap by v gene")
  
  imm_ov2.j <- repOverlap(x$data, .method = "morisita", .col="j")
  plot2.j<- vis(imm_ov2.j, .text.size = 2)+
    my_xaxis+
    sc+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample)) +ggtitle("Morisita overlap by j gene")
  
  imm_ov2.vaa <- repOverlap(x$data, .method = "morisita", .col="aa+v")
  plot2.vaa<- vis(imm_ov2.vaa, .text.size = 2)+
    my_xaxis+
    sc+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
                     limits = rev(x$meta$Sample)) +ggtitle("Morisita overlap by v+aa")
  
  #### note to self 3/10/2023
  #### may want to obtain summary statistics, e.g. mean from upper triangle. this is code to do that
  #### imm_ov4[upper.tri(imm_ov4)]
  
  #### counts > 1
  x_highCounts = x
  x_highCounts$data <- lapply(1:length(x$data), function(i) {
    filter(x$data[[i]], Clones > 1)
  }) %>% `names<-`(names(x$data))
  imm_ov3 <- repOverlap(x_highCounts$data, .method="morisita", .col="aa") 
  
  plot5<-  vis(imm_ov3, .text.size = 4, .signif.digits = 1)+
    ggtitle("Clones > 1")+
    my_xaxis+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)], 
                     limits = rev(x$meta$Sample))
  
  #### counts == 1
  x_lowCounts = x
  x_lowCounts$data <- lapply(1:length(x$data), function(i) {
    filter(x$data[[i]], Clones == 1)
  })%>% `names<-`(names(x$data))
  imm_ov4 <- repOverlap(x_lowCounts$data, .method="morisita", .col="aa")
  
  plot6<- vis(imm_ov4, .text.size = 4, .signif.digits = 1)+
    ggtitle("Clones = 1")+
    my_xaxis+
    scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)], 
                     limits = rev(x$meta$Sample))
  # Apply different analysis algorithms to the matrix of public clonotypes:
  # "mds" - Multi-dimensional Scaling
  
  # this one produces an error message
  # repOverlapAnalysis(imm_ov1, "mds")
  
  # "tsne" - t-Stochastic Neighbor Embedding -- note attempted to run this by sample but returned error "perplexity is too large for the number of samples"
  
  # plot7<- repOverlapAnalysis(imm_ov1, "tsne") %>%  vis()+theme(legend.position="NULL")
  
  ## this returns an error message
  # repOverlapAnalysis(imm_ov1, "mds+kmeans") %>% vis()
  
  pdf(paste0(output.path,"repOv.pdf"), width=12, height=8)
  print(plot1)
  print(plot2)
  print(list(plot2.nt,plot2.v,plot2.j,plot2.vaa))
  print(plot5)
  print(plot6)
  # print(plot7)
  dev.off()
  
  #this code was using clonal family data, which we've since removed
  # 
  # if (tcrbcr == "IGH") {
  #   xy = x
  #   xy$data <- xy$data %>% lapply(., function(x) x %>%   #immunarch repOverlap function won't accept the new column clone_id, so to trick it I remove the -cdr3.aa column and replace it with the clone_id column  
  #                                   select(-CDR3.aa) %>% rename(CDR3.aa = clone_id))
  #   
  #   imm_ov2.clonalfam <- repOverlap(xy$data, .method = "morisita", .col="aa")
  #   plot2.clonalfam <- vis(imm_ov2.clonalfam, .text.size = 2)+
  #     my_xaxis+
  #     sc+
  #     scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
  #                      limits = rev(x$meta$Sample)) +ggtitle("Morisita overlap by clonal family")
  #   
  #   imm_ov2.clonalfam.ds <- repOverlap(xy$data, .method = "morisita", .col="aa",.downsample=T)
  #   plot2.clonalfam.ds <- vis(imm_ov2.clonalfam.ds, .text.size = 2)+
  #     my_xaxis+
  #     sc+
  #     scale_y_discrete(labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
  #                      limits = rev(x$meta$Sample)) +ggtitle("Downsampled Morisita overlap by clonal family")
  #   
  #   pdf(paste0(output.path,"repOv_clonalfam.pdf"), width=12, height=8)
  #   print(plot2.clonalfam)
  #   print(plot2.clonalfam.ds)
  #   dev.off()
  #   
  # }
  
  #heatmaps
  if (all(unlist(imm_ov2) %in% c(0, NA))) {
    print(paste0("repertoire overlap heatmap not created for ",unique(x$meta$Patient.ID), " due to entire overlap matrix being 0"))
  } else {
  plot3 <- vis(imm_ov2, "heatmap2")
  rownames(imm_ov2) <- x$meta$tls_pub_id[match(rownames(imm_ov2), x$meta$Sample)]
  colnames(imm_ov2) <- x$meta$tls_pub_id[match(colnames(imm_ov2), x$meta$Sample)]
  plot4<- vis(imm_ov2, "heatmap2")
  
  pdf(paste0(output.path, "repOv_pHeatmap.pdf"), width=12, height=8)
  grid::grid.newpage()
  grid::grid.draw(plot3$gtable)
  grid::grid.newpage()
  grid::grid.draw(plot4$gtable)
  dev.off()
  }
}

### Gene usage
gene.usage.single.patient <- function(x) {
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  
  plot1<- geneUsage(x$data) %>% vis()+ggprism::theme_prism()+theme(axis.text.x = element_text(angle=90, size=6, hjust=1))
  plot2<- lapply(1:length(meta.subset), function (i) 
    geneUsage(x$data) %>% 
      vis(., .by = meta.subset[i], .meta = x$meta)) 
  
  pdf(paste0(output.path,"geneUsage.pdf"), width=9, height=7)
  print(plot1)
  print(plot2)
  dev.off()
  
}

### Diversity estimations
#visualize the 4 diversity functions across all samples
diversity.single.patient <- function(x){
  #set x axis parameters
  my_xaxis = list(
    theme_prism(),
    theme(axis.text.x = element_text(angle=90, vjust=0.25)),
    scale_x_discrete(
      labels= function(y) x$meta$tls_pub_id[match(y, x$meta$Sample)],
      limits = x$meta$Sample) #rev added here to accomodate coord_flip
  )
  
  #make plots
  div_chao <- repDiversity(x$data, "chao1")
  plot1<- vis(div_chao)+my_xaxis+theme(legend.position="NULL")
  # plot2<-vis(div_chao, .by="TLS_type", .meta=x$meta)+ggprism::theme_prism() ### commented out because otherwise won't run
  
  div_hill <- repDiversity(x$data, "hill")
  plot3<-vis(div_hill)
  # plot4<-vis(div_hill, .by="TLS_type", .meta=x$meta)+ggprism::theme_prism()
  
  div_d50 <- repDiversity(x$data, "d50")
  plot5<-vis(div_d50)+my_xaxis+theme(legend.position="NULL")
  # plot6<-vis(div_d50, .by="TLS_type", .meta=x$meta)+ggprism::theme_prism()
  
  div_div <- repDiversity(x$data, "div")
  plot7<-vis(div_div)+my_xaxis+theme(legend.position="NULL")
  # plot8<-vis(div_div, .by="TLS_type", .meta=x$meta)+ggprism::theme_prism()
  
  div_div <- repDiversity(coding(
    x$data)
    , "div", .q = 1, .col ="aa"
  ) #.q = 1 gives you exp shannon
  plot9<-vis(div_div) +my_xaxis+theme(legend.position="NULL")
  # plot10<-vis(div_div, .by="TLS_type", .meta=x$meta)+ggprism::theme_prism()
  
  pdf(paste0(output.path,"diversity.pdf"), width=9, height=7)
  print(plot1)
  # print(plot2)
  print(plot3)
  # print(plot4)
  print(plot5)
  # print(plot6)
  print(plot7)
  # print(plot8)
  print(plot9)
  # print(plot10)
  dev.off()
}
