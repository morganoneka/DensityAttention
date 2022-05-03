library(ggplot2)

# input_dir = "/Users/morganoneka/Box/My Stuff/CLAM/attention/attention_csv/"
# 
# for (file in list.files(input_dir))

file_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

# these are the "master files" that have all the X/Y coordinates 
cp_data <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/Chronic Pancreatitis.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
pdac_data <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/PDAC.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

for (r in 1:nrow(file_info)){
  
  if (is.na(file_info[r,"Diagnosis"] )){
    next
  }
  if (file_info[r,"Diagnosis"] != "PDAC" & file_info[r,"Diagnosis"] != "Chronic Pancreatitis"){
    next 
  }
  
  pt_number = strsplit(file_info[r,"slide_id"], "patient")[[1]][2]
  # patient_files = list.files(paste("/Users/morganoneka/Box/Phenotype Spreadsheets for Each Dx/", file_info[r,"Diagnosis"], sep=""), paste("^", pt_number, "_.*",sep=""))
  
  cell_info = data.frame()
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
  files <- list.files("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/gwr_output/", pattern="*.csv", full.names=TRUE)
  files_2 = files[grepl(paste("_", pt_number, "_", sep=""), files)]
  
  # read in ALL GWR output
  all_gwr_output <- lapply(files_2, read.table, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)
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
  ggsave(paste("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/", pt_number, ".png", sep=""), device="png")
  
  # write.table(all_in_square, paste("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/", pt_number, ".csv", sep=""), quote=FALSE, row.names=FALSE, sep=",")
  
  
}

# attention_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/1001.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
# cell_info <- read.table("/Users/morganoneka/Box/Phenotype Spreadsheets for Each Dx/Chronic Pancreatitis/1001_2_cell_seg_data.txt", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

# ggplot() + geom_point(data = attention_info, aes(x=X, y=Y, color=Attention))
# ggplot() + geom_point(data = cell_info, aes(x=Cell.X.Position, y=Cell.Y.Position))

# cp <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/Chronic Pancreatitis.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)


# allcp <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/Chronic Pancreatitis.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
# just1001 <- allcp[which(allcp$Patient == "1001"),]
# ggplot() + geom_point(data=just1001, aes(x=Cell.X.Position, y=Cell.Y.Position))

# we gotta rotate herrrrr
# attention_info$X_2 = max(attention_info$Y) - attention_info$Y
# attention_info$Y_2 = attention_info$X

# ggplot() + geom_point(data=attention_info, aes(x=X_2, y=Y_2, color=Attention)) 
# 
# xy_1 <- read.table("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/gwr_output/Chronic Pancreatitis_1001_CTL.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)
# ggplot() + geom_point(data = xy_1, aes(x=Y_min, y=X_min, color=Density))



# files <- list.files("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/gwr_output/", pattern="*.csv", full.names=TRUE)
# files_2 = files[grepl("_1001_", files)]
# # read in ALL GWR output
# all_gwr_output <- lapply(files_2, read.table, header=TRUE, stringsAsFactors=FALSE, sep=",", fill=TRUE)
# all_gwr_output_2 <- do.call(rbind, all_gwr_output)
# # swap x and y values
# all_gwr_output_2$X = all_gwr_output_2$Y_min
# all_gwr_output_2$Y = all_gwr_output_2$X_min
# 
# # get minx and y
# gwr_min_x = min(all_gwr_output_2$X)
# gwr_min_y = min(all_gwr_output_2$Y)
# 
# atn_min_x = min(attention_info$X_2)
# atn_min_y = min(attention_info$Y_2)
# 
# trans_x = gwr_min_x - atn_min_x
# trans_y = gwr_min_y - atn_min_y
# 
# # apply translation to attention score 
# attention_info$X_3 = attention_info$X_2 + trans_x
# attention_info$Y_3 = attention_info$Y_2 + trans_y
# 
# 
# # SCALING!!!
# gwr_max_x = max(all_gwr_output_2$X)
# gwr_max_y = max(all_gwr_output_2$Y)
# 
# atn_max_x = max(attention_info$X_3)
# atn_max_y = max(attention_info$Y_3)
# 
# # we need to add some logic in case there's empty space on the top or the right too
# # i think for that we can compare the gwr to the og image??? 
# 
# og_max_x = max(cell_info$Cell.X.Position)
# og_max_y = max(cell_info$Cell.Y.Position)
# 
# diff_x = og_max_x - gwr_max_x
# diff_y = og_max_y - gwr_max_y
# 
# 
# scale_x = (gwr_max_x + diff_x) / atn_max_x 
# scale_y = (gwr_max_y + diff_y) / atn_max_y
# attention_info$X_4 = attention_info$X_3 * scale_x
# attention_info$Y_4 = attention_info$Y_3 * scale_y
# 
# ggplot() + geom_tile(data = attention_info, aes(x=X_4, y=Y_4, fill=Attention)) + geom_point(data = cell_info, aes(x=Cell.X.Position, y=Cell.Y.Position))
# 
# summary(attention_info$Attention)
# 
# 
# # TODO: get regions for attention squares
# sorted_x = sort(unique(attention_info$X_4))
# x_square_size = sorted_x[2] - sorted_x[1]
# 
# sorted_y = sort(unique(attention_info$Y_4))
# y_square_size = sorted_y[2] - sorted_y[1]
# 
# # TODO: identify which points are in which squares
# 
# all_in_square = as.data.frame(matrix(nrow=0, ncol=191), stringsAsFactors=FALSE)
# for(i in 1:nrow(attention_info)){
#   if(attention_info[i, "Attention"] > 0){
#     row = attention_info[i,]
#     x = row$X_4
#     y = row$Y_4
#     
#     in_x = which(cell_info$Cell.X.Position >= x & cell_info$Cell.X.Position <= (x+x_square_size))
#     in_y = which(cell_info$Cell.Y.Position >= y & cell_info$Cell.Y.Position <= (y+y_square_size))
#     
#     in_square = cell_info[intersect(in_x,in_y),]
#     all_in_square = rbind(all_in_square, in_square)
#     
#   }
# }
# 
# ggplot() + geom_tile(data = attention_info, aes(x=X_4, y=Y_4, fill=Attention)) + geom_point(data = all_in_square, aes(x=Cell.X.Position, y=Cell.Y.Position))
# 
