
###so to create stacked barplot of clonotypes
#column 1 = clonotype (specie)
#column 2 = condition (t or b cell subset)
#column 3 = clones (value, i.e. number of clones of that category)

#First define function for stacked barplot of seurat.TLS or seurat.TLS.excludeSingletons object
stackedBP <- function(sPlusT_obj) {
  ggplot(sPlusT_obj, aes(x=reorder(TCRB_or_IGH, -countSum),
                         y=count, 
                         fill=active.ident)) + 
    geom_bar(stat="identity")+ #, colour="#000000", linewidth=0.001, show.legend = T) +
    theme_classic()+
    xlab(paste0(ifelse(repertoire=="T", "TCRB", "IGH"), " ", 
                "CDR3.aa"))+
    theme(axis.text.x = element_blank(), #element_text(size=rel(0.6), angle=90, vjust=1, hjust=1, face="bold"),
          axis.ticks.x=element_blank(),
          # legend.justification=c(1,1),
          # legend.position="right"#,
          #          legend.background = element_rect(fill = "white", color = "black")
    )+
    scale_x_discrete()+
    scale_y_continuous(expand = c(0, 0))+
    labs(fill="Cell type")+
    guides(fill=guide_legend(ncol=1))+
    scale_fill_manual(values=cluster_colors)
}

#create function prep.forBarplot, which takes an object like seurat.TLS and collapses the seurat object metadata into unique clonotypes (rather than unique cells)

# this was original code, preceded the function prepforBarplot
# s.plus.T.forBarplot <- seurat.TLS %>% 
#   select(TCRB_or_IGH, predicted.celltype.l2) %>% 
#   add_count(TCRB_or_IGH, name = 'countSum') %>% 
#   group_by(TCRB_or_IGH, predicted.celltype.l2) %>% 
#   mutate(count=n(), .before = "countSum")  %>% 
#   distinct(.keep_all=T) 

prepforBarplot <- function(x) {
  x %>% #seurat.TLS was original object) 
    select(TCRB_or_IGH, active.ident) %>% 
    add_count(TCRB_or_IGH, name = 'countSum') %>% 
    group_by(TCRB_or_IGH, active.ident) %>% 
    mutate(count=n(), .before = "countSum")  %>% 
    distinct(.keep_all=T) 
}
