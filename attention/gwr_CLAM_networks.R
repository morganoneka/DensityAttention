library(ggplot2)
library(GGMselect)
library(stringr)
library(tidyverse)
library(huge)
library(igraph)

file_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

cp_data <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/Chronic Pancreatitis.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
pdac_data <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/PDAC.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

all_gwr_in_square <- as.data.frame(matrix(nrow=0, ncol=8), stringsAsFactors=FALSE)

#TODO: want to split by PDAC and CP
# iterate over all the files
for (r in 1:nrow(file_info)){
  
  print(r)
  
  if (is.na(file_info[r,"Diagnosis"] )){
    next
  }
  if (file_info[r,"Diagnosis"] != "PDAC" & file_info[r,"Diagnosis"] != "Chronic Pancreatitis"){
    next 
  }
  
  pt_number = strsplit(file_info[r,"slide_id"], "patient")[[1]][2]
  # patient_files = list.files(paste("/Users/morganoneka/Box/Phenotype Spreadsheets for Each Dx/", file_info[r,"Diagnosis"], sep=""), paste("^", pt_number, "_.*",sep=""))
  
  # cell_info = data.frame()
  if (file_info[r,"Diagnosis"] == "Chronic Pancreatitis"){
    cell_info = cp_data[which(cp_data$Patient == pt_number),]
  } else{
    cell_info = pdac_data[which(pdac_data$Patient == pt_number),]
  }
  
  # i think we just want the last one
  if (length(unique(cell_info$Sample.Name)) > 1){
    cell_info = cell_info[which(cell_info$Sample.Name == tail(unique(cell_info$Sample.Name),1)),]
  }
  
  # get attention info from clam stuff
  if (!file.exists(paste("/Users/morganoneka/Box/My Stuff/CLAM/attention/attention_csv/patient", pt_number, ".csv", sep=""))){
    next
  }
  attention_info <- read.table(paste("/Users/morganoneka/Box/My Stuff/CLAM/attention/attention_csv/patient", pt_number, ".csv", sep=""), sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
  
  ggplot() + geom_point(data = attention_info, aes(x=X, y=Y, color=Attention))
  ggplot() + geom_point(data = cell_info, aes(x=Cell.X.Position, y=Cell.Y.Position))
  
  # rotating the attention info
  attention_info$X_2 = max(attention_info$Y) - attention_info$Y
  attention_info$Y_2 = attention_info$X
  
  ggplot() + geom_point(data=attention_info, aes(x=X_2, y=Y_2, color=Attention)) 
  
  # getting all the gwr output
  files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/gwr_output/", pattern="*.csv", full.names=TRUE)
  files_2 = files[grepl(paste("_", pt_number, "_", sep=""), files)]
  
  # read in ALL GWR output
  all_gwr_output <- lapply(files_2, function(x){
    tmp = read.table(x, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)
    tmp$fname = x
    return(tmp)
  })
  all_gwr_output_2 <- do.call(rbind, all_gwr_output)
  
  # swap x and y values
  all_gwr_output_2$X = all_gwr_output_2$Y_min
  all_gwr_output_2$Y = all_gwr_output_2$X_min
  
  # get minx and y
  gwr_min_x = min(all_gwr_output_2$X)
  gwr_min_y = min(all_gwr_output_2$Y)
  
  atn_min_x = min(attention_info$X_2)
  atn_min_y = min(attention_info$Y_2)
  
  trans_x = gwr_min_x - atn_min_x
  trans_y = gwr_min_y - atn_min_y
  
  # apply translation to attention score 
  attention_info$X_3 = attention_info$X_2 + trans_x
  attention_info$Y_3 = attention_info$Y_2 + trans_y
  
  # SCALING!!!
  gwr_max_x = max(all_gwr_output_2$X)
  gwr_max_y = max(all_gwr_output_2$Y)
  
  atn_max_x = max(attention_info$X_3)
  atn_max_y = max(attention_info$Y_3)
  
  # we need to add some logic in case there's empty space on the top or the right too
  # i think for that we can compare the gwr to the og image??? 
  
  og_max_x = max(cell_info$Cell.X.Position)
  og_max_y = max(cell_info$Cell.Y.Position)
  
  diff_x = og_max_x - gwr_max_x
  diff_y = og_max_y - gwr_max_y
  
  
  scale_x = (gwr_max_x + diff_x) / atn_max_x 
  scale_y = (gwr_max_y + diff_y) / atn_max_y
  attention_info$X_4 = attention_info$X_3 * scale_x
  attention_info$Y_4 = attention_info$Y_3 * scale_y
  
  ggplot() + geom_tile(data = attention_info, aes(x=X_4, y=Y_4, fill=Attention)) + geom_point(data = cell_info, aes(x=Cell.X.Position, y=Cell.Y.Position))
  
  # WHEW ok let's identify all the cells in the squares of interest now
  sorted_x = sort(unique(attention_info$X_4))
  x_square_size = sorted_x[2] - sorted_x[1]
  
  sorted_y = sort(unique(attention_info$Y_4))
  y_square_size = sorted_y[2] - sorted_y[1]
  
  # identify which points are in which squares
  all_in_square = as.data.frame(matrix(nrow=0, ncol=191), stringsAsFactors=FALSE)
  for(i in 1:nrow(attention_info)){
    if(attention_info[i, "Attention"] > 0){
      row = attention_info[i,]
      x = row$X_4
      y = row$Y_4
      
      in_x = which(cell_info$Cell.X.Position >= x & cell_info$Cell.X.Position <= (x+x_square_size))
      in_y = which(cell_info$Cell.Y.Position >= y & cell_info$Cell.Y.Position <= (y+y_square_size))
      
      in_square = cell_info[intersect(in_x,in_y),]
      all_in_square = rbind(all_in_square, in_square)
      
    }
  }
  
  if (nrow(all_in_square) == 0){
    next
  }
  
  ggplot() + geom_tile(data = attention_info, aes(x=X_4, y=Y_4, fill=Attention)) + geom_point(data = all_in_square, aes(x=Cell.X.Position, y=Cell.Y.Position))
  
  #TODO: for each gwr pixel, check if it's in a high attention region
  gwr_in_square = as.data.frame(matrix(nrow=0, ncol=8), stringsAsFactors=FALSE)
  for(i in 1:nrow(attention_info)){
    if(attention_info[i, "Attention"] > 0){
      row = attention_info[i,]
      x = row$X_4
      y = row$Y_4
      
      in_x = which(all_gwr_output_2$X_min >= x-10 & all_gwr_output_2$X_min <= (x+x_square_size+10))
      in_y = which(all_gwr_output_2$Y_min >= y-10 & all_gwr_output_2$Y_min <= (y+y_square_size+10))
      
      
      in_square = all_gwr_output_2[intersect(in_x,in_y),]
      gwr_in_square = rbind(gwr_in_square, in_square)
      
    }
  }
  
  all_gwr_in_square <- rbind(all_gwr_in_square, gwr_in_square)
  
  
}

# get diagnosis and cell type
relevant <- as.data.frame(str_split_fixed(all_gwr_in_square$fname, "_", 4)[,2:4], stringsAsFactors=FALSE)
relevant$V1 <- str_split_fixed(relevant$V1, "//",2)[,2]
relevant$V3 <- str_split_fixed(relevant$V3, "[.]", 2)[,1]

all_gwr_in_square <- cbind(all_gwr_in_square, relevant)
# all_gwr_in_square$key <- unlist(lapply(1:nrow(all_gwr_in_square), function(x){
#   return(paste(all_gwr_in_square[x, c("X", "Y")], collapse="_"))
# }))

no_dups <- all_gwr_in_square[!duplicated(all_gwr_in_square), c("X", "Y", "Density", "V1", "V2", "V3")]


wide_data <- no_dups %>% spread(V3, Density)
wide_data[is.na(wide_data)] <- 0


pdac_for_network <- as.matrix(wide_data[which(wide_data$V1 == "PDAC"), c("APC", "CTL","Helper T-Cell", "Other", "T-reg")])
cp_for_network <- as.matrix(wide_data[which(wide_data$V1 == "Chronic Pancreatitis"), c("APC", "CTL","Helper T-Cell", "Other", "T-reg")])

# pdac_output = selectFast(pdac_for_network)
# pdac_cor <- cor(t(pdac_for_network))

pdac_output = huge(pdac_for_network, cov.output = TRUE)
pdac_select = huge.select(pdac_output)
pdac_best_graph = as.matrix(pdac_select$refit)


cp_output = huge(cp_for_network[,c("APC", "CTL","Helper T-Cell", "Other")], cov.output = TRUE)
cp_select = huge.select(cp_output)
cp_best_graph = cp_select$refit
cp_best_graph = rbind(cp_best_graph, c(0,0,0,0))
cp_best_graph = cbind(cp_best_graph, c(0,0,0,0,0))

colnames(pdac_best_graph) = c("APC", "CTL","Helper T-Cell", "Other", "T-reg")
rownames(pdac_best_graph) = c("APC", "CTL","Helper T-Cell", "Other", "T-reg")

colnames(cp_best_graph) = c("APC", "CTL","Helper T-Cell", "Other", "T-reg")
rownames(cp_best_graph) = c("APC", "CTL","Helper T-Cell", "Other", "T-reg")

pdac_graph = graph_from_adjacency_matrix(pdac_best_graph, mode="undirected", diag=FALSE)
cp_graph = graph_from_adjacency_matrix(cp_best_graph, mode="undirected", diag=FALSE)

plot(pdac_graph)
plot(cp_graph)
