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

files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/", pattern=".csv", full.names=TRUE)
diagnosis_info <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")

PDAC_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "PDAC"), "slide_id"])
CP_patients = unlist(diagnosis_info[which(diagnosis_info$Diagnosis == "Chronic Pancreatitis"), "slide_id"])

# TODO split files vector into two vectors- one for pdac and one for cp
patientnames = unlist(lapply(strsplit(files,"//"), function(x) return( paste("patient",strsplit(x[[2]], "[.]")[[1]][1], sep="") )))

PDAC_files = files[patientnames %in% unlist(PDAC_patients)]
CP_files = files[patientnames %in% unlist(CP_patients)]

# TODO create two data framesbased on the region info
PDAC_all <- do.call(rbind, lapply(PDAC_files, function(x){
  tmp = read.table(x, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")
  if (nrow(tmp) > 0){
    tmp$file = x
    return(tmp)
  } else{
    return(NULL)
  }
  
}))

CP_all <- do.call(rbind, lapply(CP_files, function(x){
  tmp = read.table(x, header=TRUE, fill=TRUE, stringsAsFactors=FALSE, sep=",")
  if (nrow(tmp) > 0){
    tmp$file = x
    return(tmp)
  } else{
    return(NULL)
  }
  
}))

# TODO OFFSET!!! 
for (i in 1:length(CP_files)){
  idx = which(CP_all$file == CP_files[i])
  CP_all[idx, "Cell.X.Position"] = CP_all[idx, "Cell.X.Position"] + 2000*i
  CP_all[idx, "Cell.Y.Position"] = CP_all[idx, "Cell.Y.Position"] + 2000*i
}

for (i in 1:length(PDAC_files)){
  idx = which(PDAC_all$file == PDAC_files[i])
  PDAC_all[idx, "Cell.X.Position"] = PDAC_all[idx, "Cell.X.Position"] + 2000*i
  PDAC_all[idx, "Cell.Y.Position"] = PDAC_all[idx, "Cell.Y.Position"] + 2000*i
}

for (data in c(CP_all, PDAC_all)){
    
    
    data_in <- data[!duplicated(data[,c("Cell.X.Position","Cell.Y.Position")]),]
    data_in <- data_in[,c("Cell.X.Position", "Cell.Y.Position", "CellType")]
    
    data_in$expr = 1
    
    non_gene_columns <- c("Cell.X.Position", "Cell.Y.Position", "CellType")
    x_col = "Cell.X.Position"
    y_col = "Cell.Y.Position"
    
    
    print("Making Giotto object")
    # go <- createGiottoObject(transpose(data_in[,colnames(data_in)[!(colnames(data_in) %in% non_gene_columns)]]), data_in[,c(x_col,y_col)])
    go <- createGiottoObject(transpose(data_in[,c("expr","expr")]), data_in[,c(x_col,y_col)])
    print("Giotto object made")
    
    # with_annot <- addGeneMetadata(go, data.frame(gene_ID = c("CD8", "FoxP3", "PD.L1", "CD4", "Treg","APC","Epithelial", "Tcell")))
    with_network <- createSpatialNetwork(go)
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
    
    # filename = paste(patient_and_sample, gsub(" ", "_", tissue_condition), sep="_")
    

    # write.table(to_write_out[order(-abs(to_write_out$enrichm)),],file=fullfile("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/", paste(filename, "interaction_enrichment.csv", sep="_")), sep=",", row.names=FALSE, quote=FALSE)
    # write.table(to_write_out[order(-abs(to_write_out$enrichm)),],file=fullfile("/Users/morganoneka/Box/My Stuff/CLAM/attention/giotto_output/", paste(filename, "interaction_enrichment.csv", sep="_")), sep=",", row.names=FALSE, quote=FALSE)
  }
      
      
      
      

    
    











