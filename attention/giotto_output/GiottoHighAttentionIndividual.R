library("matlab")
library("optparse")
library("Giotto")
library("reshape2")
library("ggtern")
library("akima")
library("ggalluvial")
library("igraph")
library("stringr")
library(tidyr)


# for each group
output_dir="/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/individual/"
files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/", pattern=".csv", full.names=TRUE)

# patients <- unique(str_split_fixed(files, "_", 2)[,1])

# read in files from each patient
for (file in files){
  
  # if (file == "Sample_1201991_trans.csv" | file == "Sample_1220211_trans.csv" | file ==  "Sample_1229818_trans.csv"){
  #   print("Problematic file")
  #   next
  # }
  # if (file != "Sample_1201991_trans.csv" & file != "Sample_1220211_trans.csv" & file !=  "Sample_1229818_trans.csv"){
  #   print("Already ran Giotto")
  #   next
  # }
  
  # print(paste(group, patient))
  # patient_and_sample = paste(str_split_fixed(file,"_",3)[,1:2], collapse="_")
  
  # for (pfile in patient_files){
  data_in <- read.table(file, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")

  if (nrow(data_in) == 0){
    next
  }
  
  data_in <- data_in[!duplicated(data_in[,c("Cell.X.Position","Cell.Y.Position")]),]
  
  
  
  data_in$expr = 1
  
  #subset only cells where Markers are valid
  celltypes = unique(data_in$CellType)
  # celltypes_relevant = celltypes[grepl("pSTAT", celltypes)]
  
  # data_relevant = data_in[which(data_in$Markers %in% celltypes_relevant),]
  
  #Tsplit based on Type (so tumor vs. advancing edge vs. necrosis)
  # data_split = split(data_relevant, data_relevant$Type, drop=FALSE)
  # print(length(data_split))

    
    filename = strsplit(tail(strsplit(file,"/")[[1]],1), "[.]")[[1]][1]
    if(file.exists(fullfile(output_dir, paste(filename, "interaction_enrichment.csv", sep="_")))){
      print("File Exists")
      
    } else{
      #TODO: change all of this to work with current data
      # non_gene_columns <- c("Markers", "X", "Y", "Type")
      x_col = "Cell.X.Position"
      y_col = "Cell.Y.Position"
      
      
      print("Making Giotto object")
      # go <- createGiottoObject(transpose(data_in[,colnames(data_in)[!(colnames(data_in) %in% non_gene_columns)]]), data_in[,c(x_col,y_col)])
      go <- tryCatch({createGiottoObject(transpose(data_in[,c("expr","expr")]), data_in[,c(x_col,y_col)])}, error=function(cond){
        print(cond)
        return(NA)
      })
      
      if (is.na(go)){
        next
      }
      print("Giotto object made")
      
      # with_annot <- addGeneMetadata(go, data.frame(gene_ID = c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")))
      with_network <- tryCatch({createSpatialNetwork(go)}, error=function(cond){
        print(cond)
        return(NA)
      })
      
      if(is.na(with_network)){
        next
      }
      
      with_spatialgrid <- createSpatialGrid(with_network, sdimx_stepsize = 50, sdimy_stepsize = 50)
      print("Spatial Grid made")
      
      
      with_phenotype <- with_spatialgrid
      with_phenotype@cell_metadata$phenotype <- data_in$CellType
      print("Reassigned phenotype")
      
      # tryCatch(, error = next)
      cell_prox_enrich <- cellProximityEnrichment(with_phenotype, cluster_column="phenotype")
      
      # i <- 0
      # while (TRUE){
      #   if (i > 9){
      #     print("reached thresh")
      #     break
      #   }
      #   tryCatch(cell_prox_enrich <- cellProximityEnrichment(with_phenotype, cluster_column="phenotype"),
      #            error = function(e){
      #              print("error, trying again")
      #              print(paste("i is ", i))
      #               i<-i+1
      #            }
      #            )
      # }
      
      print("Proximity enrichment")
      
      # plot cell enrichment
      #TODO: make sure plotting is working
      cell_prox_enrich$enrichm$cell_1 <- str_split_fixed(cell_prox_enrich$enrichm$unified_int,"--",2)[,1]
      cell_prox_enrich$enrichm$cell_2 <- str_split_fixed(cell_prox_enrich$enrichm$unified_int,"--",2)[,2]
      
      to_write_out <- cell_prox_enrich$enrichm[,c("enrichm", "p_higher_orig", "p_lower_orig", "PI_value", "cell_1", "cell_2")]
      
     
      
      write.table(to_write_out[order(-abs(to_write_out$enrichm)),],file=fullfile(output_dir, paste(filename, "interaction_enrichment.csv", sep="_")), sep=",", row.names=FALSE, quote=FALSE)
    
    
    
    
  }
  
  
}










