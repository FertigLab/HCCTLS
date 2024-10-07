TLS_composition <- function(query){ #use if there is only one seurat object to check
  meta_subset = adaptiveTLS$meta[grep(query, adaptiveTLS$meta$Sample),]
  TLS_sampNames <- meta_subset$Sample
  list <- adaptiveTLS$data[TLS_sampNames]
  TLSrep <- bind_rows(list,.id="TLS") %>% select(TLS, Clones, Proportion, CDR3.nt, CDR3.aa)  
  # TLSrep %>% head
  TLSrep %>% group_by(TLS) %>% summarise(total_proportion=sum(Proportion)) #this shows that the Proportion column doesnt add up to 100%, because non-productive sequences have been excluded, so we have to make a new one based on Clones column
  TLSrep = TLSrep %>% group_by(TLS) %>% mutate(new_prop = Clones/sum(Clones)) %>% ungroup()
  TLSrep %>% group_by(TLS) %>% summarise(total_new_prop = sum(new_prop)) #this doublechecks that this adds up to 100 by TLS now
  TLSrep %>% head #look at new_prop, which should be slightly higher than Proportion
  TLSrep = select(TLSrep, -Proportion) %>% rename(., Proportion = new_prop) #drops old Proportion then renames new one as Proportion
  TLSrep %>% head
  #####################################################
  # Build single cell table for that same patient
  scRep <- seurat@meta.data %>% filter(Patient==query) %>% select(active.ident, TCRB_or_IGH, TCRB_or_IGH2)
  scRep <- mutate_all(scRep, funs(replace(., .=='NA', NA))) # some TCRB_or_IGH have character string "NA" rather than NA, so change these to form that is.na will recognize and omit
  
  # scRep$TCRB_or_IGH %>% is.na() %>% table
  
  #remove rows with NA in TCRB_or_IGH or TCRB_or_IGH2
  scRep.dropNA <- scRep[!with(scRep,is.na(TCRB_or_IGH)&is.na(TCRB_or_IGH2)),]
  #scRep.dropNA$TCRB_or_IGH %>% is.na() %>% table
  
  # scRep.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n()) %>% View
  
  scRep.dropNA.summarised <- scRep.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n())
  
  scRep.dropNA.summarised.wider <- scRep.dropNA.summarised %>% pivot_wider(., names_from = active.ident, values_from = TCRB_or_IGH_n)
  scRep.dropNA.summarised.wider$subsets = rowSums( !is.na(scRep.dropNA.summarised.wider[,2:ncol(scRep.dropNA.summarised.wider)]))
  # scRep.dropNA.summarised.wider %>% View
  scRep_final = scRep.dropNA.summarised.wider %>% 
    group_by(TCRB_or_IGH) %>% summarise_each((funs(sum))) %>% 
    # mutate(top_r=apply(.[,2:(ncol(.)-1)], 1, function(x) names(x)[which.max(x)]))
    mutate(top_r=apply(.[,names(.) %in% unique(scRep$active.ident)], 1, function(x) names(x)[which.max(x)]))
  
  # View(scRep_final)
  
  ####################################################
  #make df by joining the two objects above (TLS data and single cell data)
  df <- left_join(TLSrep, select(scRep_final, c(TCRB_or_IGH,top_r)), 
                  by=c("CDR3.aa" = "TCRB_or_IGH"))
  df$top_r = df$top_r %>% replace_na("unmatched")
  df$top_r <- factor(x=df$top_r, levels=c(as.character(levels.manual), "unmatched")) #set levels and add "unmatched" to levels
  # cluster_colors = c(cluster_colors,"lightgray")
  # names(cluster_colors)[length(cluster_colors)] = "unmatched"
  df$expanded = if_else(df$Clones==1, "singleton", "expanded")
  df$known = if_else(df$top_r == "unmatched", 0, 1)
  df = merge(df,meta_subset[,c("Sample","TLS_type")], by.x="TLS",by.y="Sample") 
  df$TLS = df$TLS %>% str_replace("-T_TCRB","")

  df.expanded = df[df$expanded=="expanded",] 
  df.expanded = df.expanded %>% group_by(TLS) %>% mutate(id = row_number())
  
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_1.pdf"), height=12, width=16)
  n_to_test = c(10, 15, 20, 25, 30,40,50,100)
  TLS.labs = meta_subset$tls_pub_id
  names(TLS.labs) = meta_subset$Sample
  for (n_clonotypes in n_to_test){
    ptop15 <- df.expanded %>% group_by(TLS) %>% 
      slice_head(n=n_clonotypes) %>%   
      ggplot(., aes(x=id, 
                    y=Proportion, fill=top_r))+
      geom_bar(stat="identity",color="black",size=0.25)+
      ggprism::theme_prism()+
      theme(axis.line.x = element_blank(),
            axis.text.x = element_blank(),#(angle=45),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text=element_text(face="bold")
      )+
      labs(fill="")+
      scale_fill_manual(values=cluster_colors)+
      # theme_minimal()+
      facet_wrap(TLS_type ~ TLS,nrow=2,
                 labeller = labeller(TLS = TLS.labs))+
      ggtitle(paste0("Inferred composition of Top ",n_clonotypes, " clonotypes by TLS"))
    print(ptop15)
  }
  dev.off()
  
  
  #now do the same but just look at recovery of expanded clones
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_2.pdf"), height=12, width=16)
  TLS.labs = meta_subset$tls_pub_id
  names(TLS.labs) = meta_subset$Sample
  for (i in c(1, 0.5, 0.2, 0.1,0.05,0.01)){
    data2 = df.expanded %>% slice_head(prop=i) %>% 
      group_by(TLS) %>%
      mutate(new_prop = Clones/sum(Clones)) %>%  #recalculates proportion to proportion of expanded only
      select(-Proportion) %>% rename(., Proportion = new_prop) %>%  ungroup()
    data2=data2 %>% 
      group_by(TLS,top_r,known,TLS_type) %>%
      summarise(Proportion_sum=sum(Proportion),
                Clones_sum=sum(Clones)) %>% ungroup()
    data2 = data2 %>% group_by(TLS) %>% 
      arrange(.,TLS,desc(known),desc(Proportion_sum)) %>%  #arrange in desired order
      mutate(ymax = cumsum(Proportion_sum),
             ymin = c(0,head(ymax,n=-1)))
    
    p0.1 <- ggplot(data2,aes(ymax=ymax, ymin=ymin))+
      # geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
      geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      scale_fill_manual(values=c(cluster_colors,"expanded"="#1B191999", "singleton" = "#80818099"))+
      labs(x = NULL, y = NULL) +
      theme(axis.text.x = element_blank(),#(angle=45),
            axis.ticks.x = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text=element_text(face="bold")
      )+
      facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)")) 
    print(p0.1)
    
    p0.1.1 <- ggplot(data2,aes(x=TLS,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))+
      scale_x_discrete(
        labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
        # limits = meta_subset$Sample
      )
    print(p0.1.1)
    
    p0.1.1.1 <- ggplot(data2,aes(x=TLS,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))+
      scale_x_discrete(
        labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
        # limits = meta_subset$Sample
      )+
      facet_wrap(~TLS_type, scales="free_x")
    print(p0.1.1.1)
    
    p0.1.2 <- ggplot(data2,aes(x=TLS_type,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))
    print(p0.1.2)
  }
  
  data3=df[!df$top_r %in% ("unmatched"),] %>%  #another iteration, drop the unmatched and plot all knowns, expanded or not expanded
    group_by(TLS) %>%
    mutate(new_prop = Clones/sum(Clones)) %>%  #recalculates proportion to proportion of expanded only
    select(-Proportion) %>% rename(., Proportion = new_prop) %>%  ungroup()
  data3=data3 %>% 
    group_by(TLS,top_r,expanded,TLS_type) %>%
    summarise(Proportion_sum=sum(Proportion),
              Clones_sum=sum(Clones)) %>% ungroup()
  data3 = data3 %>% group_by(TLS) %>% 
    arrange(.,TLS,expanded,
            desc(Proportion_sum)) %>%  #arrange in desired order
    mutate(ymax = cumsum(Proportion_sum),
           ymin = c(0,head(ymax,n=-1)))
  p0.2 <- ggplot(data3,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2)
  
  p0.2.1 <- ggplot(data3,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    # coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs),
               nrow=1)+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2.1)
  
  p0.2.2 <- ggplot(data3,aes(x=TLS,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))+
    scale_x_discrete(
      labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
      # limits = meta_subset$Sample
    )  
  print(p0.2.2)
  
  p0.2.2.1 <- ggplot(data3,aes(x=TLS,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))+
    scale_x_discrete(
      labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
      # limits = meta_subset$Sample
    )+
    facet_wrap(~TLS_type, scales ="free_x")
  print(p0.2.2.1)
  
  p0.2.3 <- ggplot(data3,aes(x=TLS_type,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2.3)
  
  
  #make collapsed versions of df by collapsing unique CDR3s to make for faster plotting -- I doublechecked this and it comes to the same thing
  df_collapsed <- df %>% group_by(TLS,top_r,known,expanded) %>%
    summarise(Proportion_sum=sum(Proportion),
              Clones_sum=sum(Clones))
  
  data = df_collapsed %>% group_by(TLS) %>%
    arrange(.,TLS,expanded,desc(known),desc(Proportion_sum)) %>%  #arrange in desired order
    mutate(ymax = cumsum(Proportion_sum),
           ymin = c(0,head(ymax,n=-1)))
  # alternate code for making data object, not used because becomes impossible to plot
  # data = df %>% group_by(TLS) %>%
  #   arrange(.,desc(Proportion))%>%  #arrange in desired order
  #   mutate(ymax = cumsum(Proportion),
  #          ymin = c(0,head(ymax,n=-1)))
  
  p0 <- ggplot(data,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+ggtitle("TLS TCR phenotype matching")
  print(p0)
  
  dev.off()
  
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_3.pdf"), height=12, width=16)
  
  # ##plot 1 and 2 - plots of entire repertoire, without alpha
  # p1 <- ggplot(df_collapsed, aes(x = TLS, y=Proportion_sum, fill=top_r))+
  #   geom_bar(position="stack", stat="identity")+
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
  # print(p1)
  # 
  # p2 <- ggplot(df_collapsed, aes(x = TLS, y=Clones_sum, fill=top_r))+ #collapsed 
  #   geom_bar(position="stack", stat="identity")+
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p2)
  
  # #plot 3 and 4, plots of entire repertoire, with alpha
  # p3<-ggplot(df_collapsed, aes(x = TLS, y=Proportion_sum, fill = top_r, alpha = expanded)) +
  #   geom_bar(position="stack", stat="identity")+#, color='black')+
  #   scale_alpha_manual(values=c(1,0.4)) +
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p3)  
  # 
  # p4<-ggplot(df_collapsed, aes(x = TLS, y=Clones_sum, fill = top_r, alpha = expanded)) +
  #   geom_bar(position="stack", stat="identity")+#, color='black')+
  #   scale_alpha_manual(values=c(1,0.4)) +
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p4)
  
  #plots using subsetting for top X% clones, without alpha (5,7,9,11) and with alpha (6,8,10,12) 
  for (i in c(1,0.5, 0.1, 0.05, 0.01,0.005)) {
    df_top <- df %>%
      group_by(TLS) %>%
      slice(seq(n()*i)) #0.1 = top 10%
    # df_top=df
    df_top = df_top %>% group_by(TLS,top_r,known,expanded,TLS_type) %>%
      summarise(Proportion=sum(Proportion),
                Clones=sum(Clones))
    
    my_xaxis = list(
      scale_fill_manual(values=cluster_colors),
      scale_y_continuous(expand = c(0.01,0.01)),
      ggprism::theme_prism(),
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS")),
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)),
      scale_x_discrete(
        labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
        # limits = meta_subset$Sample
      )
    )
    
    p1 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="stack", stat="identity")+
      my_xaxis
    
    print(p1)
    
    p2 <- ggplot(df_top, aes(x = TLS, y=Clones, fill=top_r))+ #collapsed 
      geom_bar(position="stack", stat="identity")+
      my_xaxis
    
    print(p2)
    
    p3 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      my_xaxis
    print(p3)
    
    p3.1 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      my_xaxis+
      facet_wrap(~TLS_type, scales="free_x")
    print(p3.1)
    
    p4 <- ggplot(df_top, aes(x = TLS_type, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=cluster_colors)+
      scale_y_continuous(expand = c(0.01,0.01))+
      ggprism::theme_prism()+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS, by TLS type"))+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
    print(p4)
    
    p5<- ggplot(df_top[df_top$top_r!="unmatched",], aes(x = TLS_type, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=cluster_colors)+
      scale_y_continuous(expand = c(0.01,0.01))+
      ggprism::theme_prism()+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS, by TLS type"))+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
    print(p5)
  }
  dev.off()
}
## Add an alpha value to a colour
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

TLS_composition_til_pbmc <- function(query, seurat_obj_TIL, seurat_obj_PBMC){
  meta_subset = adaptiveTLS$meta[grep(query, adaptiveTLS$meta$Sample),]
  TLS_sampNames <- meta_subset$Sample
  list <- adaptiveTLS$data[TLS_sampNames]
  TLSrep <- bind_rows(list,.id="TLS") %>% select(TLS, Clones, Proportion, CDR3.nt, CDR3.aa)  
  # TLSrep %>% head
  TLSrep %>% group_by(TLS) %>% summarise(total_proportion=sum(Proportion)) #this shows that the Proportion column doesnt add up to 100%, because non-productive sequences have been excluded, so we have to make a new one based on Clones column
  TLSrep = TLSrep %>% group_by(TLS) %>% mutate(new_prop = Clones/sum(Clones)) %>% ungroup()
  TLSrep %>% group_by(TLS) %>% summarise(total_new_prop = sum(new_prop)) #this doublechecks that this adds up to 100 by TLS now
  TLSrep %>% head #look at new_prop, which should be slightly higher than Proportion
  TLSrep = select(TLSrep, -Proportion) %>% rename(., Proportion = new_prop) #drops old Proportion then renames new one as Proportion
  TLSrep %>% head
  #####################################################
  # Build single cell table for that same patient
  scRep <- seurat_obj_TIL@meta.data %>% filter(Patient==query) %>% select(active.ident, TCRB_or_IGH, TCRB_or_IGH2)
  scRep <- mutate_all(scRep, funs(replace(., .=='NA', NA))) # some TCRB_or_IGH have character string "NA" rather than NA, so change these to form that is.na will recognize and omit
  
  # scRep$TCRB_or_IGH %>% is.na() %>% table
  
  #remove rows with NA in TCRB_or_IGH or TCRB_or_IGH2
  scRep.dropNA <- scRep[!with(scRep,is.na(TCRB_or_IGH)&is.na(TCRB_or_IGH2)),]
  # scRep.dropNA$TCRB_or_IGH %>% is.na() %>% table
  
  # scRep.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n()) %>% View
  
  scRep.dropNA.summarised <- scRep.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n())
  
  scRep.dropNA.summarised.wider <- scRep.dropNA.summarised %>% pivot_wider(., names_from = active.ident, values_from = TCRB_or_IGH_n)
  # scRep.dropNA.summarised.wider$subsets = rowSums( !is.na(scRep.dropNA.summarised.wider[,2:ncol(scRep.dropNA.summarised.wider)]))
  scRep.dropNA.summarised.wider$subsets = rowSums( !is.na(scRep.dropNA.summarised.wider[,names(scRep.dropNA.summarised.wider) %in% unique(scRep$active.ident)]))
  # scRep.dropNA.summarised.wider %>% View
  scRep_final = scRep.dropNA.summarised.wider %>% 
    group_by(TCRB_or_IGH) %>% summarise_each((funs(sum))) %>% 
    # mutate(top_r=apply(.[,2:(ncol(.)-1)], 1, function(x) names(x)[which.max(x)]))
    mutate(top_r=apply(.[,names(.) %in% unique(scRep$active.ident)], 1, function(x) names(x)[which.max(x)]))
  
  
  #add to top_r column prefix that identifies that this came from TIL
  # scRep_final$top_r = paste0("TIL_",scRep_final$top_r) 
  scRep_final$top_r = paste0(scRep_final$top_r,"_TIL") 
  # View(scRep_final)
  
  ####################################################
  # Build a second single cell table for that same patient, now using the second object seurat_obj_PBMC
  scRep2 <- seurat_obj_PBMC@meta.data %>% filter(Patient==query) %>% select(active.ident, TCRB_or_IGH, TCRB_or_IGH2)
  scRep2 <- mutate_all(scRep2, funs(replace(., .=='NA', NA))) # some TCRB_or_IGH have character string "NA" rather than NA, so change these to form that is.na will recognize and omit
  
  # scRep2$TCRB_or_IGH %>% is.na() %>% table
  
  #remove rows with NA in TCRB_or_IGH or TCRB_or_IGH2
  scRep2.dropNA <- scRep2[!with(scRep2,is.na(TCRB_or_IGH)&is.na(TCRB_or_IGH2)),]
  # scRep2.dropNA$TCRB_or_IGH %>% is.na() %>% table
  
  # scRep2.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n()) %>% View
  
  scRep2.dropNA.summarised <- scRep2.dropNA %>%  group_by(TCRB_or_IGH, active.ident) %>% summarise(TCRB_or_IGH_n = n())
  
  scRep2.dropNA.summarised.wider <- scRep2.dropNA.summarised %>% pivot_wider(., names_from = active.ident, values_from = TCRB_or_IGH_n)
  # scRep2.dropNA.summarised.wider$subsets = rowSums( !is.na(scRep2.dropNA.summarised.wider[,2:ncol(scRep2.dropNA.summarised.wider)]))
  scRep2.dropNA.summarised.wider$subsets = rowSums(!is.na(scRep2.dropNA.summarised.wider[,names(scRep2.dropNA.summarised.wider) %in% unique(scRep2$active.ident)]))
  
  # scRep2.dropNA.summarised.wider %>% View
  scRep2_final = scRep2.dropNA.summarised.wider %>% 
    group_by(TCRB_or_IGH) %>% summarise_each((funs(sum))) %>% 
    # mutate(top_r=apply(.[,2:(ncol(.)-1)], 1, function(x) names(x)[which.max(x)]))
    mutate(top_r=apply(.[,names(.) %in% unique(scRep$active.ident)], 1, function(x) names(x)[which.max(x)]))
  
  #add to top_r column prefix that identifies that this came from PBMC
  # scRep2_final$top_r = paste0("PBMC_",scRep2_final$top_r) 
  scRep2_final$top_r = paste0(scRep2_final$top_r,"_PBMC")
  
  # View(scRep2_final)
  
  
  ####################################################
  #make df by joining the two objects above (TLS data and single cell data)
  df_TLS_TIL <- left_join(TLSrep, select(scRep_final, c(TCRB_or_IGH,top_r)), 
                  by=c("CDR3.aa" = "TCRB_or_IGH")) %>% 
    rename(top_r_TIL = top_r)
  df <- left_join(df_TLS_TIL, select(scRep2_final, c(TCRB_or_IGH,top_r)), 
                  by=c("CDR3.aa" = "TCRB_or_IGH")) %>% 
    rename(top_r_PBMC = top_r)
  # check recovery
  # df$top_r_TIL %>% is.na %>% table  #recovers 980 clonotypes across all TLS (duplicates if same clonotype in multiple tls)
  # df$top_r_PBMC %>% is.na %>% table #recovers 2031 clonotypes across all TLS

  # df[is.na(df$top_r_TIL)&is.na(df$top_r_PBMC),] %>% nrow#this is the df where we didn't match anything
  # df[!is.na(df$top_r_TIL)&is.na(df$top_r_PBMC),] %>% nrow #this is the df where we matched TIL but not PBMC
  # df[is.na(df$top_r_TIL)&!is.na(df$top_r_PBMC),] %>% View #this is the df where we matched PBMC but not TIL  
  # df[!is.na(df$top_r_TIL)&!is.na(df$top_r_PBMC),] %>% View #this is the df where we matched PBMC but not TIL  
  
  #next we need to make a final top_r column, by preferentially taking TIL phenotype first, and if that isn't available take PBMC phenotype
  df$top_r = 
    if_else(!is.na(df$top_r_TIL)==T, 
            df$top_r_TIL,
            if_else(!is.na(df$top_r_PBMC)==T,
                    df$top_r_PBMC,
                    "unmatched"))
  #alternate code
    # paste0(df$top_r_TIL,"_", df$top_r_PBMC)
    #then replace "NA_NA" with "unmatched"
    # df[df$top_r=="NA_NA",]$top_r <- "unmatched"
    #and drop the _NA or NA_ from the columns where there was a partial match
    # df$top_r %>% is.na %>% table #this pre-post check makes sure I didn't mess anything up with the str_replace in next line
    # df$top_r = str_replace(df$top_r, "NA_|_NA","")
  # df$top_r %>% table
  # df$top_r %>% table %>% sum
  df$top_r = str_replace(df$top_r, "_PBMC|_TIL","")
  # df$top_r %>% table
  # df$top_r %>% table %>% sum
  
  df$top_r <- factor(x=df$top_r, levels=c(as.character(levels.manual), "unmatched"))
  ####
  # this code I previously used to create separate cluster oclors for both TIL and PBMC, commenting out now since I collapse these two groups ##
  # if i want to revert to this, i would comment out the line above, with the factor assignment for top_r and and assign levels with preservation of the _TIL and _PBMC suffixes
  ####
  # levels_top_r=sort(unique(df$top_r))
  # df$top_r = factor(df$top_r, levels=levels_top_r)
  
  # cluster_colors_orig = cluster_colors
  # cluster_colors_TIL = cluster_colors %>% `names<-`(., paste0(names(.),"_TIL"))
  # cluster_colors_PBMC = add.alpha(cluster_colors,0.5) %>% `names<-`(., paste0(names(.),"_PBMC"))#since i'm using the TIL and PBMC suffix, I have to make a new cluster_colors object that uses those names instead
  # new_cc = c()
  # for (i in 1:length(cluster_colors_TIL)) {
  #   new_cc = c(new_cc, cluster_colors_TIL[i], cluster_colors_PBMC[i])
  # }
  # cluster_colors=c(cluster_colors, "lightgray")
  # names(cluster_colors)[length(cluster_colors)] <- "unmatched" #name last item, "lightgray", as unmatched
  # levels(cluster_colors) = cluster_colors
  # cluster_colors
  #####
 

  df$expanded = if_else(df$Clones==1, "singleton", "expanded")
  df$known = if_else(df$top_r == "unmatched", 0, 1)
  df = merge(df,meta_subset[,c("Sample","TLS_type")], by.x="TLS",by.y="Sample") 
  df$TLS = df$TLS %>% str_replace("-T_TCRB","")
  
  df.expanded = df[df$expanded=="expanded",] 
  df.expanded = df.expanded %>% group_by(TLS) %>% mutate(id = row_number())
  
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_1_TIL_and_PBMC.pdf"), height=5, width=16)
  n_to_test = c(10, 15, 20, 25, 30,40,50,100)
  TLS.labs = meta_subset$tls_pub_id
  names(TLS.labs) = meta_subset$Sample
  for (n_clonotypes in n_to_test){
  ptop15 <- df.expanded %>% group_by(TLS) %>% 
    slice_head(n=n_clonotypes) %>%   
    ggplot(., aes(x=id, 
                  y=Proportion, fill=top_r))+
    geom_bar(stat="identity",color="black",size=0.25)+
    ggprism::theme_prism()+
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
          )+
    labs(fill="")+
    scale_fill_manual(values=cluster_colors)+
    # theme_minimal()+
    facet_wrap(TLS_type ~ TLS,nrow=1,
               labeller = labeller(TLS = TLS.labs))+
    ggtitle(paste0("Inferred phenotype of Top ",n_clonotypes, " clonotypes by TLS"))+
    xlab("Top TCR\u03b2 clonotypes")
  print(ptop15)
  }
  dev.off()
  

  #now do the same but just look at recovery of expanded clones
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_2_TIL_and_PBMC.pdf"), height=12, width=16)
  TLS.labs = meta_subset$tls_pub_id
  names(TLS.labs) = meta_subset$Sample
  for (i in c(1, 0.5, 0.2, 0.1,0.05,0.01)){
    data2 = df.expanded %>% slice_head(prop=i) %>% 
      group_by(TLS) %>%
      mutate(new_prop = Clones/sum(Clones)) %>%  #recalculates proportion to proportion of expanded only
      select(-Proportion) %>% rename(., Proportion = new_prop) %>%  ungroup()
    data2=data2 %>% 
      group_by(TLS,top_r,known,TLS_type) %>%
      summarise(Proportion_sum=sum(Proportion),
                Clones_sum=sum(Clones)) %>% ungroup()
    data2 = data2 %>% group_by(TLS) %>% 
      arrange(.,TLS,desc(known),desc(Proportion_sum)) %>%  #arrange in desired order
      mutate(ymax = cumsum(Proportion_sum),
             ymin = c(0,head(ymax,n=-1)))
    
    p0.1 <- ggplot(data2,aes(ymax=ymax, ymin=ymin))+
      # geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
      geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      scale_fill_manual(values=c(cluster_colors,"expanded"="#1B191999", "singleton" = "#80818099"))+
      labs(x = NULL, y = NULL) +
      theme(axis.text.x = element_blank(),#(angle=45),
            axis.ticks.x = element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            strip.background = element_blank(),
            strip.text=element_text(face="bold")
      )+
      facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)")) 
    print(p0.1)
    
    p0.1.1 <- ggplot(data2,aes(x=TLS,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))+
      scale_x_discrete(
        labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
        # limits = meta_subset$Sample
      )
    print(p0.1.1)
    
    p0.1.1.1 <- ggplot(data2,aes(x=TLS,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))+
      scale_x_discrete(
        labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
        # limits = meta_subset$Sample
      )+
      facet_wrap(~TLS_type, scales="free_x")
    print(p0.1.1.1)
    
    p0.1.2 <- ggplot(data2,aes(x=TLS_type,y=Proportion_sum, fill=top_r))+
      # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
      # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
      geom_bar(stat="identity",position="fill")+
      scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
      ))+
      labs(y = "Proportion") +
      ggprism::theme_prism()+
      theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of TCR\u03B2 in TLS (expanded only)"))
    print(p0.1.2)
  }
  
  data3=df[!df$top_r %in% ("unmatched"),] %>%  #another iteration, drop the unmatched and plot all knowns, expanded or not expanded
    group_by(TLS) %>%
    mutate(new_prop = Clones/sum(Clones)) %>%  #recalculates proportion to proportion of expanded only
    select(-Proportion) %>% rename(., Proportion = new_prop) %>%  ungroup()
  data3=data3 %>% 
    group_by(TLS,top_r,expanded,TLS_type) %>%
    summarise(Proportion_sum=sum(Proportion),
              Clones_sum=sum(Clones)) %>% ungroup()
  data3 = data3 %>% group_by(TLS) %>% 
    arrange(.,TLS,expanded,
            desc(Proportion_sum)) %>%  #arrange in desired order
    mutate(ymax = cumsum(Proportion_sum),
           ymin = c(0,head(ymax,n=-1)))
  p0.2 <- ggplot(data3,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2)
  
  p0.2.1 <- ggplot(data3,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    # coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs),
               nrow=1)+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2.1)
  
  p0.2.2 <- ggplot(data3,aes(x=TLS,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))+
    scale_x_discrete(
      labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
      # limits = meta_subset$Sample
    )  
  print(p0.2.2)
  
  p0.2.2.1 <- ggplot(data3,aes(x=TLS,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))+
    scale_x_discrete(
      labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
      # limits = meta_subset$Sample
    )+
    facet_wrap(~TLS_type, scales ="free_x")
  print(p0.2.2.1)
  
  p0.2.3 <- ggplot(data3,aes(x=TLS_type,y=Proportion_sum, fill=top_r))+
    # geom_rect(aes(fill=top_r, xmin=1, xmax=2))+
    # geom_rect(aes(fill=top_r, xmin=2, xmax=3))+coord_polar(theta="y")+xlim(c(0,3))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=c(cluster_colors#,"expanded"="#1B191999", "singleton" = "#80818099"
    ))+
    labs(y = "Proportion") +
    ggprism::theme_prism()+
    theme(axis.text.x = element_text(angle=45, hjust=1,vjust=1))+
    ggtitle(paste0("Phenotype of all matched TCR\u03B2 in TLS (unmatched excluded)"))
  print(p0.2.3)
  
  
  #make collapsed versions of df by collapsing unique CDR3s to make for faster plotting -- I doublechecked this and it comes to the same thing
  df_collapsed <- df %>% group_by(TLS,top_r,known,expanded) %>%
    summarise(Proportion_sum=sum(Proportion),
              Clones_sum=sum(Clones))
  
  data = df_collapsed %>% group_by(TLS) %>%
    arrange(.,TLS,expanded,desc(known),desc(Proportion_sum)) %>%  #arrange in desired order
    mutate(ymax = cumsum(Proportion_sum),
           ymin = c(0,head(ymax,n=-1)))
  # alternate code for making data object, not used because becomes impossible to plot
  # data = df %>% group_by(TLS) %>%
  #   arrange(.,desc(Proportion))%>%  #arrange in desired order
  #   mutate(ymax = cumsum(Proportion),
  #          ymin = c(0,head(ymax,n=-1)))
  
  p0 <- ggplot(data,aes(ymax=ymax, ymin=ymin))+
    geom_rect(aes(fill=expanded, xmin=1, xmax=2))+
    scale_fill_manual(values=c("expanded"="#1B191999", "singleton" = "#80818099"))+
    new_scale_fill()+
    geom_rect(aes(fill=top_r, xmin=2, xmax=3))+
    scale_fill_manual(values=c(cluster_colors))+
    coord_polar(theta="y")+xlim(c(0,3))+
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),#(angle=45),
          axis.ticks.x = element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          strip.text=element_text(face="bold")
    )+
    facet_wrap(~TLS,labeller = labeller(TLS = TLS.labs))+ggtitle("TLS TCR phenotype matching")
  print(p0)
  
  dev.off()
  
  pdf(paste0(output.path,"TLS_composition_plot_", query, "_3_TIL_and_PBMC.pdf"), height=12, width=16)
  
  # ##plot 1 and 2 - plots of entire repertoire, without alpha
  # p1 <- ggplot(df_collapsed, aes(x = TLS, y=Proportion_sum, fill=top_r))+
  #   geom_bar(position="stack", stat="identity")+
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
  # print(p1)
  # 
  # p2 <- ggplot(df_collapsed, aes(x = TLS, y=Clones_sum, fill=top_r))+ #collapsed 
  #   geom_bar(position="stack", stat="identity")+
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p2)
  
  # #plot 3 and 4, plots of entire repertoire, with alpha
  # p3<-ggplot(df_collapsed, aes(x = TLS, y=Proportion_sum, fill = top_r, alpha = expanded)) +
  #   geom_bar(position="stack", stat="identity")+#, color='black')+
  #   scale_alpha_manual(values=c(1,0.4)) +
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p3)  
  # 
  # p4<-ggplot(df_collapsed, aes(x = TLS, y=Clones_sum, fill = top_r, alpha = expanded)) +
  #   geom_bar(position="stack", stat="identity")+#, color='black')+
  #   scale_alpha_manual(values=c(1,0.4)) +
  #   scale_fill_manual(values=cluster_colors)+
  #   scale_y_continuous(expand = c(0.01,0.01))+
  #   ggprism::theme_prism()+
  #   ggtitle("TLS T cell phenotypes recovered by cross-referencing single cell data, by TLS")+theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
  # print(p4)
  
  #plots using subsetting for top X% clones, without alpha (5,7,9,11) and with alpha (6,8,10,12) 
  for (i in c(1,0.5, 0.1, 0.05, 0.01,0.005)) {
        df_top <- df %>%
      group_by(TLS) %>%
      slice(seq(n()*i)) #0.1 = top 10%
        # df_top=df
    df_top = df_top %>% group_by(TLS,top_r,known,expanded,TLS_type) %>%
    summarise(Proportion=sum(Proportion),
              Clones=sum(Clones))
    
    my_xaxis = list(
      scale_fill_manual(values=cluster_colors),
        scale_y_continuous(expand = c(0.01,0.01)),
        ggprism::theme_prism(),
        ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS")),
        theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)),
        scale_x_discrete(
          labels= function(y) meta_subset$tls_pub_id[match(y, meta_subset$Sample)]#,
          # limits = meta_subset$Sample
          )
    )
    
    p1 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="stack", stat="identity")+
      my_xaxis
      
    print(p1)
    
    p2 <- ggplot(df_top, aes(x = TLS, y=Clones, fill=top_r))+ #collapsed 
      geom_bar(position="stack", stat="identity")+
      my_xaxis
    
    print(p2)
    
    p3 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      my_xaxis
    print(p3)
    
    p3.1 <- ggplot(df_top, aes(x = TLS, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
     my_xaxis+
      facet_wrap(~TLS_type, scales="free_x")
    print(p3.1)
    
    p4 <- ggplot(df_top, aes(x = TLS_type, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=cluster_colors)+
      scale_y_continuous(expand = c(0.01,0.01))+
      ggprism::theme_prism()+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS, by TLS type"))+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
    print(p4)
    
    p5<- ggplot(df_top[df_top$top_r!="unmatched",], aes(x = TLS_type, y=Proportion, fill=top_r))+
      geom_bar(position="fill", stat="identity")+
      scale_fill_manual(values=cluster_colors)+
      scale_y_continuous(expand = c(0.01,0.01))+
      ggprism::theme_prism()+
      ggtitle(paste0("Inferred phenotype of top ", i*100, "% of all TCR\u03B2 in TLS, by TLS type"))+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1)) 
    print(p5)
  }
  dev.off()
  
  
#final plot for all unmatched  
  df_drop_unmatched <- df %>% filter(top_r != "unmatched") %>% 
    group_by(TLS_type,top_r) %>%
    summarise(Proportion=sum(Proportion),
              Clones=sum(Clones))
  
  p<-ggplot(df_drop_unmatched, aes(x = TLS_type, y=Proportion, fill=top_r))+
    geom_bar(position="fill", stat="identity",color="white")+
    scale_fill_manual(values=cluster_colors)+
    scale_y_continuous(expand = c(0.01,0.01),labels=scales::percent_format(suffix=""))+
    ggprism::theme_prism()+
    ggtitle(paste0("Inferred phenotype of TCR\u03B2 \nin TLS (Patient ", query, ")"))+
      theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
            title=element_text(hjust=0.5))+
    xlab("TLS type")+ylab("Proportion (%)")+
    guides(fill=guide_legend(ncol=1))

  pdf(paste0(output.path,"TLS_composition_plot_", query, "_4_for_pub_TIL_and_PBMC.pdf"),width=5)
  print(p)
  dev.off()
  
  # df_drop_unmatched_collapsed <- df_drop_unmatched %>%
  #   mutate(top_r_collapsed = str_replace(top_r, "_PBMC|_TIL","")) %>%
  #   group_by(TLS_type,top_r_collapsed) %>%
  #   summarise(Proportion=sum(Proportion),
  #             Clones=sum(Clones))
  # df_drop_unmatched_collapsed$top_r_collapsed <- factor(df_drop_unmatched_collapsed$top_r_collapsed, levels = levels.manual)
  # 
  # p<-ggplot(df_drop_unmatched_collapsed, aes(x = TLS_type, y=Proportion, fill=top_r_collapsed))+
  #   geom_bar(position="fill", stat="identity",color="white")+
  #   scale_fill_manual(values=cluster_colors_orig)+
  #   scale_y_continuous(expand = c(0.01,0.01),labels=scales::percent_format(suffix=""))+
  #   ggprism::theme_prism()+
  #   ggtitle(paste0("Inferred phenotype of TCR\u03B2 \nin TLS of OT6 (unmatched excluded)"))+
  #   theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1),
  #         title=element_text(hjust=0.5))+
  #   xlab("TLS type")+ylab("Proportion (%)")+
  #   guides(fill=guide_legend(ncol=1))
  # 
  # pdf(paste0(output.path,"TLS_composition_plot_", query, "_5_for_pub.pdf"),width=5)
  # print(p)
  # dev.off()
  
}
