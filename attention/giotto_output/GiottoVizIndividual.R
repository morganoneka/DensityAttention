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
library(tidyverse)


id_enriched_depleted <- function(x){
  x$InxType = "Neither"
  x$InxType[x$p_higher_orig <= 0.05 ]  ="Enrichment"
  x$InxType[x$p_lower_orig <= 0.05 ]  ="Depletion"
  # x$cell_1 <- unlist(lapply(x$cell_1,only_positives))
  # x$cell_2 <- unlist(lapply(x$cell_2,only_positives))
 return(x)
} 


files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/individual/", pattern="*.csv", full.names=TRUE)
diagnosis_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")


PDAC_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "PDAC"), "slide_id"])
CP_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "Chronic Pancreatitis"), "slide_id"])

# TODO split files vector into two vectors- one for pdac and one for cp
patientnames = unlist(lapply(list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/individual/", pattern="*.csv", full.names=FALSE),
                      function(x){
                        paste("patient", strsplit(x, "_")[[1]][1], sep="")
                      }))

PDAC_files = files[patientnames %in% unlist(PDAC_patients)]
CP_files = files[patientnames %in% unlist(CP_patients)]

# TODO create two data framesbased on the region info
PDAC_all <- lapply(PDAC_files, function(x){
  tmp = read.table(x, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")
  if (nrow(tmp) > 0){
    tmp$file = x
    return(id_enriched_depleted(tmp))
  } else{
    return(NULL)
  }

})

CP_all <- lapply(CP_files, function(x){
  tmp = read.table(x, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")
  if (nrow(tmp) > 0){
    tmp$file = x
    return(id_enriched_depleted(tmp))
  } else{
    return(NULL)
  }

})

all_data <- lapply(c(PDAC_files, CP_files), function(x){
    tmp = read.table(x, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")
    if (nrow(tmp) > 0){
      tmp$file = x
      return(tmp)
    } else{
      return(NULL)
    }

  })

all_data_clean <- lapply(all_data, id_enriched_depleted)


# consensus_edge
diff_network_object <- create_diff_network_object(all_data, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")
CP_object <- create_diff_network_object(CP_all, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")
PDAC_object <- create_diff_network_object(PDAC_all, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")

consensus_edge(CP_object, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")
consensus_edge(PDAC_object, node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")


make_piecharts <- function(x){
  node1="cell_1"
  node2="cell_2"
  diff_network_object = create_diff_network_object(all_data_clean[x], node1="cell_1", node2="cell_2", edge_weight="enrichm", edge_color="InxType")
  node_xy_df <- diff_network_object[[2]]
  
  edge_list_labels = c(node1,node2)
  
  edge_list <- do.call(rbind,lapply(diff_network_object[[1]], FUN = function(x){x[,c(edge_list_labels, "InxType")]} ))
  
  edge_list_in_order <- edge_list[which(edge_list[,node1] <= edge_list[,node2]),]
  edge_list_out_of_order <- edge_list[which(edge_list[,node1] > edge_list[,node2]),]
  colnames(edge_list_out_of_order) <- c(node2,node1,  "InxType")
  
  
  occurrences <- melt(table(rbind(edge_list_in_order, edge_list_out_of_order)))
  # occurrences <- occurrences[which(occurrences$value >0),]
  
  # occurrences <- merge(occurrences, node_xy_df, by.x=node1, by.y="Node", all.x=TRUE, all.y=FALSE)
  # occurrences <- merge(occurrences, node_xy_df, by.x=node2, by.y="Node", all.x=TRUE, all.y=FALSE)
  
  celltypes <- sort(unlist(diff_network_object[[2]]$Node))
  coords <- as.data.frame(cbind(celltypes, 1:length(celltypes)), stringsAsFactors = FALSE)
  coords$V2 = as.numeric(coords$V2)
  
  combo = merge(occurrences, coords, by.x="cell_1", "celltypes", all=FALSE)
  combo = merge(combo, coords, by.x="cell_2", "celltypes", all=FALSE)
  
  colnames(combo) <- c("cell_1", "cell_2", "InxType", "value", "X", "Y")
  # combo[is.na(combo$value),"value"] = 0
  
  combo_flipped  = merge(occurrences, coords, by.x="cell_2", "celltypes", all=FALSE)
  combo_flipped = merge(combo_flipped, coords, by.x="cell_1", "celltypes", all=FALSE)
  
  colnames(combo_flipped) <- c("cell_1", "cell_2", "InxType", "value", "X", "Y")
  # combo_flipped[is.na(combo_flipped$value),"value"] = 0
  
  combo = rbind(combo, combo_flipped)
  
  combo = combo[which(combo$value >0),]
  
  num_patients = length(x)
  
  combo$cell_1 = as.character(combo$cell_1)
  combo$cell_2 = as.character(combo$cell_2)
  
  for (cell1 in celltypes){
    for (cell2 in celltypes){
      for (inxtype in c("Enrichment", "Depletion")){
        
        w = which(combo$cell_1 == cell1 & combo$cell_2 == cell2 & combo$InxType == inxtype)
        # print(w)
        if (length(w) == 0){
          combo = rbind(combo, c(cell1,cell2,inxtype,0,coords[which(coords$celltypes == cell1), "V2"],coords[which(coords$celltypes == cell2), "V2"]), stringsAsFactors=FALSE)
        }
      }
    }
  }
  
  colnames(combo) <- c("cell_1", "cell_2", "InxType", "value", "X", "Y")
  
  combo$X = as.numeric(combo$X)
  combo$Y = as.numeric(combo$Y)
  combo$InxType = as.character(combo$InxType)
  combo$value = as.numeric(combo$value)
  combo$value = combo$value / num_patients
  
  combo_clean <- dcast(unique(combo), cell_1 + cell_2 ~ InxType, value.var="value")
  combo_clean = merge(combo_clean, coords, by.x="cell_1", "celltypes", all=FALSE)
  combo_clean = merge(combo_clean, coords, by.x="cell_2", "celltypes", all=FALSE)
  colnames(combo_clean) <- c("cell_1", "cell_2", "Depletion", "Enrichment", "Neither", "X", "Y")
  combo_clean$Neither = 1 - combo_clean$Depletion - combo_clean$Enrichment
  combo_clean$radius = 0.25
  
  
  # print(combo_clean)
  
  ggplot() + geom_scatterpie(aes(x=X, y=Y, r=radius), data=combo_clean, cols=c("Depletion", "Enrichment", "Neither")) + coord_equal() +
    scale_x_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) +
    scale_y_discrete(limits = as.character(0:(length(celltypes)-1) + 0.5), breaks=as.character(0:(length(celltypes)-1) + 0.5), labels=celltypes) + 
    theme_pubr() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size=9), axis.text.y = element_text(size=9), axis.title.x = element_blank(), axis.title.y = element_blank())  + 
    scale_fill_manual(values=c("#3e7ced", "#ed523e", "#aaaaaa")) + labs(fill = "Interaction Type")
  
  return(combo_clean)
  
}

# PDAC
PDAC_vals = make_piecharts(1:length(PDAC_files))

# CP
CP_vals = make_piecharts((length(PDAC_files)+1):length(files))



test_output = as.data.frame(do.call(rbind, lapply(1:nrow(PDAC_vals), function(x){
  unlist(lapply(c("Depletion", "Enrichment", "Neither"), function(inxtype){
    prop.test(c(PDAC_vals[x,inxtype]*length(PDAC_files), CP_vals[x,inxtype]*length(CP_files)), c(length(PDAC_files), length(CP_files)))$p.value
  }))
})))

colnames(test_output)<-c("Depletion", "Enrichment", "Neither")
test_output$Interaction = unlist(lapply(1:nrow(PDAC_vals), function(x) {
    return(paste(PDAC_vals[x, c("cell_1", "cell_2")], collapse="_"))
  }))

test_output[which(test_output$Depletion <= 0.05),]
test_output[which(test_output$Enrichment <= 0.05),]
test_output[which(test_output$Neither <= 0.05),]

test_output_exact_test = as.data.frame(do.call(rbind, lapply(1:nrow(PDAC_vals), function(x){
  return(c(paste(PDAC_vals[x, c("cell_1", "cell_2")], collapse="_"), fisher.test(rbind(PDAC_vals[x,c("Depletion", "Enrichment", "Neither")]*length(PDAC_files), CP_vals[x,c("Depletion", "Enrichment", "Neither")]*length(CP_files)))$p.value))
})), stringsAsFactors = FALSE)

colnames(test_output_exact_test) <- c("Interaction", "Pval")
test_output_exact_test[which(test_output_exact_test$Pval <= 0.05),]

# test_output_exact_test$Pval <- as.numeric(test_output_exact_test$Pval)

# prop.test(PDAC_vals[1,c("Depletion", "Enrichment", "Neither")], CP_vals[1,c("Depletion", "Enrichment", "Neither")])

# ggarrange(make_piecharts(indices[[1]]), make_piecharts(indices[[2]]),make_piecharts(indices[[3]]), make_piecharts(indices[[4]]),make_piecharts(indices[[5]]), common.legend = TRUE, legend = "right", labels=c("NASH Advanced", "NASH Minimal", "HCV Advanced", "HCV Minimal", "HCV Control"))
# ggsave("/Users/morganoneka/Box/Nash Data/Results/Giotto/FinalFigs/PieCharts.pdf", device="pdf", width=12, height=8)


region_dir = "/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/"


# Getting the treg/treg enriched 
idx_treg = lapply(1:length(PDAC_all), function(x){
  idx = which(PDAC_all[[x]][,"cell_1"] == "T-reg" & PDAC_all[[x]][,"cell_2"] == "T-reg")
  return(idx != 0 & PDAC_all[[x]][idx, "InxType"] == "Enrichment")
})

PDAC_files[which(idx_treg == TRUE)]

pt1019 = read.table("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/1019.csv", sep=",", fill=TRUE, header=TRUE, stringsAsFactors = FALSE)
pt420 = read.table("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/420.csv", sep=",", fill=TRUE, header=TRUE, stringsAsFactors = FALSE)

pt64 = read.table("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/64.csv", sep=",", fill=TRUE, header=TRUE, stringsAsFactors = FALSE)


ggplot(pt1019, aes(x=Cell.X.Position, y=Cell.Y.Position, color=CellType))  + geom_point()
ggplot(pt420, aes(x=Cell.X.Position, y=Cell.Y.Position, color=CellType))  + geom_point()

ggplot(pt64, aes(x=Cell.X.Position, y=Cell.Y.Position, color=CellType))  + geom_point()

idx_tregtumor = lapply(1:length(PDAC_all), function(x){
  idx = which(PDAC_all[[x]][,"cell_1"] == "T-reg" & PDAC_all[[x]][,"cell_2"] == "Tumor")
  return(idx != 0 & PDAC_all[[x]][idx, "InxType"] == "Depletion")
})
PDAC_files[which(idx_tregtumor == TRUE)]






