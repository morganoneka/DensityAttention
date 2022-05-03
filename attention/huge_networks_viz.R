library("matlab")
library("optparse")
library("Giotto")
library("reshape2")
library(ggplot2)
library("ggtern")
library("akima")
library("ggalluvial")
library("igraph")
library("stringr")
library(tidyr)
library(ggforce)
library(reshape2)
library(ggrepel)
library(ggnewscale)
library(scatterpie)
library(ggpubr)

files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/huge_networks/", pattern=".txt", full.names=TRUE)
diagnosis_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")

PDAC_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "PDAC"), "slide_id"])
CP_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "Chronic Pancreatitis"), "slide_id"])

# TODO split files vector into two vectors- one for pdac and one for cp
patientnames = unlist(lapply(strsplit(files,"//"), function(x) return( paste("patient",strsplit(x[[2]], "[.]")[[1]][1], sep="") )))

PDAC_files = files[patientnames %in% unlist(PDAC_patients)]
CP_files = files[patientnames %in% unlist(CP_patients)]


edge_lists <- lapply(c(PDAC_files, CP_files), function(fname){
  # fname = "/Users/morganoneka/Box/My Stuff/CLAM/attention/huge_networks/19.txt"
  
  con <- file(fname,"r")
  col_names <- readLines(con,n=1)
  close(con)
  col_names = str_replace(col_names, "Helper T-Cell", "Helper-T-Cell")
  
  
  adj_mx <- read.table(fname, header=FALSE, fill=TRUE, stringsAsFactors=FALSE, skip=1)
  colnames(adj_mx) = unlist(strsplit(col_names, " "))
  g <- graph.adjacency(as.matrix(adj_mx))
  edge_list = get.edgelist(g)
  mx_to_use = cbind(edge_list, rep(1, nrow(edge_list)))
  colnames(mx_to_use) <- c("cell_1", "cell_2", "enrich")
  
  return(as.data.frame(mx_to_use, stringsAsFactors = FALSE))
})

diffnet_PDAC = create_diff_network_object(edge_lists[1:length(PDAC_files)], node1="cell_1", node2="cell_2", edge_weight="enrichm")

consensus_edge(diffnet_PDAC, node1="cell_1", node2="cell_2", edge_weight="enrichm")

diffnet_CP = create_diff_network_object(edge_lists[(length(PDAC_files)+1):(length(PDAC_files)+length(CP_files))], node1="cell_1", node2="cell_2", edge_weight="enrichm")

consensus_edge(diffnet_CP, node1="cell_1", node2="cell_2", edge_weight="enrichm")



