library(ggplot2)
library(GGMselect)
library(stringr)
library(tidyverse)
library(huge)
library(igraph)

file_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", sep=",", header=TRUE, fill=TRUE, stringsAsFactors=FALSE)

data_dir = "/Users/morganoneka/Box/My Stuff/CLAM/attention/gwr_output/"
output_dir = "/Users/morganoneka/Box/My Stuff/CLAM/attention/huge_networks/"
file_list <- list.files(data_dir, pattern="*.csv", full.names = FALSE)

all_celltypes =  c("APC", "CTL","Helper T-Cell", "Other", "T-reg")

# iterate over list
for (patient in unique(file_info$case_id)){
  patient_number = str_split_fixed(patient, "patient", 2)[,2]
  
  if (file.exists(paste(output_dir, patient_number, ".txt", sep=""))){
    print("file exists")
    next
  }
  
  patient_files = file_list[grepl(paste(".*_", patient_number, "_.*", sep=""), file_list)]
  
  if(length(patient_files) == 0){
    next
  }
  
  data_in = do.call(cbind, lapply(patient_files, function(x){
    tmp = read.table(paste(data_dir,x,sep=""), header=TRUE, fill=TRUE, sep=",", stringsAsFactors = FALSE)
    celltype = str_split_fixed(str_split_fixed(x, "_", 3)[,3], "[.]", 2)[,1]
    colnames(tmp) <- c("X_min", "Y_min", celltype, "X_max", "Y_max")
    # return(as.data.frame(tmp[,c("X_min", "Y_min", celltype)]))
    return(tmp[,celltype])
  }))
  
  if (ncol(data_in) <= 1){
    next
  }
  
  cell_types_patient = lapply(patient_files, function(x){
    tmp = read.table(paste(data_dir,x,sep=""), header=TRUE, fill=TRUE, sep=",", stringsAsFactors = FALSE)
    return(str_split_fixed(str_split_fixed(x, "_", 3)[,3], "[.]", 2)[,1])
  })
  
  colnames(data_in) <- cell_types_patient
  # excluded_celltypes = setdiff(all_celltypes, colnames(data_in))
  # 
  # data_in = cbind(data_in, matrix(0, nrow=nrow(data_in), ncol=length(excluded_celltypes)))
  # colnames(data_in) <- c(cell_types_patient, excluded_celltypes)
  # 
  output_network =  huge(data_in, cov.output = TRUE)
  output_best = as.matrix(huge.select(output_network)$refit)
  
  colnames(output_best) = cell_types_patient
  rownames(output_best) = cell_types_patient
  
  # TODO: when writing out, check if patient name is in PDAC or CP
  write.table(output_best, paste(output_dir, patient_number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, quote=FALSE)
}
