library(stringr)

groups = c("Chronic Pancreatitis", "IPMN", "MCN", "PanIN", "PDAC")

group_patient_df = data.frame(nrow=0, ncol=2)
colnames(group_patient_df) <- c("File", "Diagnosis")

for (group in groups){
  output_dir=paste("/Users/morganoneka/Box/My Stuff/GraphClassification/GiottoOutput/", group, sep="")
  
  files <- list.files(path=output_dir, pattern="txt")
  diagnosis_df <- as.data.frame(cbind(files, rep(group, length(files))))
  colnames(diagnosis_df) <- c("File", "Diagnosis")
  
  group_patient_df <- rbind(group_patient_df, diagnosis_df)
}



group_patient_df$Patient <- str_split_fixed(group_patient_df$File, "_", 2)[,1]

group_patient_df$HasPatches = FALSE
group_patient_df[which(group_patient_df$Patient %in% unlist(files_with_patches$Patient)), "HasPatches"] = TRUE

#####


files_with_patches <- read.table("/Users/morganoneka/Box/My Stuff/CLAM/all_patients.txt", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
files_with_patches$Patient <- str_split_fixed(str_split_fixed(files_with_patches$V1,".jpg",2)[,1], "patient",2)[,2]
files_with_patches$case_id <- unlist(lapply(files_with_patches$Patient, function(x){
  return(paste("patient",x,sep=""))
}))

files_with_patches$slide_id <- files_with_patches$case_id 


all_files <- list.files("/Users/morganoneka/Documents/Grad School/Code/NeighborhoodCoordination-master/PancreaticData/gwr_output")
all_files_split <- as.data.frame(str_split_fixed(all_files,"_",3)[,1:2])
colnames(all_files_split) <- c("Diagnosis", "Patient")

all_files_trim <- all_files_split[!duplicated(all_files_split$Patient),]

merged <- merge(files_with_patches, all_files_trim, by="Patient", all.x=TRUE, all.y=FALSE)
merged_subset <- merged[which(merged$Diagnosis == "Chronic Pancreatitis" | merged$Diagnosis == "PDAC"),]

write.table(merged[,c("case_id", "slide_id", "Diagnosis")], "/Users/morganoneka/Box/My Stuff/CLAM/info_CP_PDAC.csv", row.names=FALSE, quote=FALSE, sep=",")


