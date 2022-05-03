source("/Users/morganoneka/Documents/PersonalProjects/DifferentialNetworkVis/DifferentialNetworkVis.R")
library(ggpubr)
library(stringr)
library(igraph)
library(ggplot2)
library(ggforce)
library(reshape2)
library(ggrepel)
library(ggnewscale)
library(scatterpie)

# adapted from my dif network viz code
plot_patient_alt <-function(diff_network_object,idx=1, node1="E1", node2="E2", edge_weight="Weight", edge_color="Type"){
  
  edge_info <- diff_network_object[[1]][[idx]]
  node_xy_df <- diff_network_object[[2]]
  
  # get node1 coords
  edge_info <- merge(edge_info, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  edge_info <- merge(edge_info, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  
  names(edge_info)[names(edge_info) == edge_weight] <- 'Weight'
  names(edge_info)[names(edge_info) == edge_color] <- 'Type'
  
  
  edge_info$size <- 1.5
  
  
  # split based on edge_type
  split_data <- split(edge_info,edge_info$Type)
  g <- ggplot()
  
  # TODO: generate colors based on different groups - need high and low value
  colors <- data.frame(high=c("red", "gray50"), low=c("blue", "gray90"))
  for (i in 1:length(split_data)){
    # logic for self loops
    splits <- split_loops(split_data[[i]])
    
    edges_loops <- splits[[1]]
    edges_nonloops <- splits[[2]]
    
    # g = g 
    # + scale_color_gradient(high="#ff0000", low="#0000ff")
    # scale_color_gradient(high=colors[i,"high"],low=colors[i,"low"])
    
    # TODO: need to generalize 
    # thanks to: https://github.com/eliocamp/ggnewscale 
    if (i == 1){
      print("i = 1")
      edges_loops[,"Depleted Weight"] <- edges_loops$Weight
      edges_nonloops[,"Depleted Weight"] <- edges_nonloops$Weight
      g = g + geom_curve(
        aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Depleted Weight`, size=size), alpha=0.8,
        data = edges_nonloops) +
        geom_curve(
          aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Depleted Weight`, size=size), alpha=0.8,  curvature = 50,angle = 270,
          data = edges_loops) + scale_color_gradient(low="#191970",high="#ADD8E6", guide = guide_legend(order = 1)) 
    } else{
      print("i = 2")
      edges_loops[,"Enriched Weight"] <- edges_loops$Weight
      edges_nonloops[,"Enriched Weight"] <- edges_nonloops$Weight
      g = g + new_scale_color() + geom_curve(
        aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Enriched Weight`, size=size), alpha=0.8,
        data = edges_nonloops) +
        geom_curve(
          aes(x = X.x, y = Y.x, xend = X.y, yend = Y.y, color=`Enriched Weight`, size=size), alpha=0.8,  curvature = 50,angle = 270,
          data = edges_loops)   + scale_color_gradient(high="#9A2A2A",low="#FAA0A0", guide = guide_legend(order = 2)) 
    }
  }
  
  
  g = g + geom_circle(aes(x0 = X, y0 = Y, r = radius, fill = Node), data=node_xy_df) + 
    geom_label_repel(aes(x=X,y=Y,label=Node),hjust=0, vjust=0, force=15, size=3, data=node_xy_df) +
    guides(fill=guide_legend(title="Cell Type")) + guides(size = FALSE) +
    ylim(c(0, max(node_xy_df$Y+10))) + xlim(c(0, max(node_xy_df$X+10))) +
    coord_fixed() + theme_void()
  
  g
}

id_enriched_depleted <- function(x){
  x$InxType = "Neither"
  x$InxType[x$p_higher_orig <= 0.05 ]  ="Enrichment"
  x$InxType[x$p_lower_orig <= 0.05 ]  ="Depletion"
  # x$cell_1 <- unlist(lapply(x$cell_1,only_positives))
  # x$cell_2 <- unlist(lapply(x$cell_2,only_positives))
  x = x[which(x$InxType != "Neither"),]
  # relevant_cell_types <- c("CD16+CD163+MAC387+", "CD163+MAC387+", "CD68+", "CD14+CD68+CD163+MAC387+", "CD14+CD68+CD163+", "MAC387+", "Negative")
  # relevant_cell_types <- c("CD3+ (480) pSTAT+", "CD3+ (480) pSTAT-", "CD163+ (690) pSTAT+", "CD163+ (690) pSTAT-", "CD68+ CD163+ (570 690) pSTAT+", 
  #                          "CD68+ CD163+ (570 690) pSTAT-", "CD11c+ CD68+ CD163+ (780 570 690) pSTAT-", "CD11c+ CD68+ CD163+ (780 570 690) pSTAT+", 
  #                          "CD11c+ CD68+ CD163+ (780 570 690) pSTAT+", "CD11c+ CD68+ CD163+ (780 570 690) pSTAT-", "CD11c+ CD163+ (780 690) pSTAT-", 
  #                          "CD68+ (570) pSTAT+", "CD68+ (570) pSTAT-")
  # x = x[which(x$cell_1 %in% relevant_cell_types & x$cell_2 %in% relevant_cell_types),]
  # x$cell_1 = gsub("( )+", " ", gsub("\\([0-9 ]*\\)", "", unlist(x$cell_1)))
  # x$cell_2 = gsub("( )+", " ", gsub("\\([0-9 ]*\\)", "", unlist(x$cell_2)))
  
  return(x)
} 

fix_celltype_split <- function(x){
  problematic_rows = which(startsWith(x$cell_2, "-"))
  
  for (i in problematic_rows){
    x[i,"cell_1"] = paste(x[i,"cell_1"], "-", sep="")
    x[i,"cell_2"] = str_split_fixed(x[i,"cell_2"], "-", 2)[2]
  }
  
  return(x)
}

file_list = c("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/CP_interaction_enrichment.csv", "/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/PDAC_interaction_enrichment.csv")
dfs <- lapply(unlist(file_list), read.table, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep="," )
dfs_x <- lapply(dfs, id_enriched_depleted)

# dfs_x_fixed <- lapply(dfs_x, fix_celltype_split)

objects <- create_diff_network_object(dfs_x, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")

plot_patient_alt(objects, 1, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")  + theme(plot.title = element_text(size=22)) + theme(plot.subtitle = element_text(size=18))
ggsave("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/CP.png", device="png", height=8, width=10)
plot_patient_alt(objects, 2, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")  + theme(plot.title = element_text(size=22)) + theme(plot.subtitle = element_text(size=18))
ggsave("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/PDAC.png", device="png", height=8, width=10)

# patient_names <- unlist(lapply(unlist(file_list), function(x){
#   return(str_split_fixed(tail(strsplit(x, "/")[[1]],1), "inter",2)[,1]      )
# }))
# 
# patient_df <- as.data.frame(str_split_fixed(patient_names, "_", 3), stringsAsFactors=FALSE)[,2:3]
# colnames(patient_df) <- c("Sample", "Region")
# patient_df$Region = gsub("_", " ", unlist(patient_df$Region))
# 
# for (sample_number in unique(patient_df$Sample)){
#   pt_idx = c()
#   
#   for (idx in which(patient_df$Sample == sample_number)){
#     if (nrow(objects[[1]][idx][[1]]) > 0){
#       pt_idx = c(pt_idx, idx)
#     }
#   }
#   
#   
#   plots_pt_1 <- lapply(pt_idx, function(x){
#     return (plot_patient_alt(objects, x, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType"))  + theme(plot.title = element_text(size=22)) + theme(plot.subtitle = element_text(size=18))
#   })
#   
#   fig1 <- ggarrange(plotlist=plots_pt_1,common.legend = TRUE, legend = "right", labels=patient_df[pt_idx, "Region"]) 
#   annotate_figure(fig1, top=text_grob(paste("Cell Pair Interaction Profiles:", sample_number), size=16)) + theme(plot.margin = unit(c(1,1,1,1), "cm"))
#   
#   ggsave(paste("/Users/morganoneka/Box/Spatial analysis data/Giotto Analysis/NetworkViz/", sample_number, ".png", sep=""), device="png", height=8, width=10)
# }
# 
# 
# 
