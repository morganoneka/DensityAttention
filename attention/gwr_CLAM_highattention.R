library(spatstat)
library(sp)
library(reshape2)
library(GWmodel)
library(ggplot2)

thres = 1e-4
# kde_grid_size = c(25,25)
square_size=10
mincount=5

# the cell types we will compare to tumor
nonref_celltypes = c("APC", "CTL", "Helper T-Cell", "Other", "T-reg")

# where we're saving all this out to
output_dir = "/Users/morganoneka/Box/My Stuff/CLAM/attention/gwr_output/"

# the files to iterate over
files <- list.files("/Users/morganoneka/Box/My Stuff/CLAM/attention/high_attention_regions/", pattern=".csv", full.names=TRUE)

for (file in files){
  print(file)
  # read in data
  diagnosis_data <- read.table(file, header=TRUE, sep=",", fill=TRUE, stringsAsFactors=FALSE)
  
  # get each individual patient
  split_data <- split(diagnosis_data, diagnosis_data$Sample.Name)
  group_name <- diagnosis_data[1, "Group"]
  
  # iterate over each patient
  for (analysis_data in split_data){
    
    patient_id = analysis_data[1,"Patient"]
    print(patient_id)
    
    # look at each cell type
    for (nonrefcell in nonref_celltypes){
      print(nonrefcell)
      if (file.exists(paste(output_dir, paste(group_name, patient_id, nonrefcell, sep="_"), ".csv", sep=""))){
        print("file already exists")
        next()
      }
      # table(analysis_data$CellType)
      analysis_data$CellOfInterest = ""
      analysis_data[which(analysis_data$CellType == "Tumor"), "CellOfInterest"] = "Ref"
      analysis_data[which(analysis_data$CellType == nonrefcell), "CellOfInterest"] = "Other"
      
      if (sum(unlist(analysis_data$CellOfInterest) == "Ref")<mincount | 
          sum(unlist(analysis_data$CellOfInterest) == "Other")<mincount){
        next
      }
      
      dataPoints=data.frame(Cell.X.Position=(analysis_data$Cell.X.Position),
                            Cell.Y.Position=(analysis_data$Cell.Y.Position), 
                            CellOfInterest=analysis_data$CellOfInterest)
      
      ref = dataPoints[dataPoints$CellOfInterest=='Ref',]
      nonref = dataPoints[dataPoints$CellOfInterest=='Other',]
      
      dataPoints = rbind(ref,nonref)
      
      xmin = min(analysis_data$Cell.X.Position)
      xmax = max(analysis_data$Cell.X.Position)
      ymin = min(analysis_data$Cell.Y.Position)
      ymax = max(analysis_data$Cell.Y.Position)
      
      xstep = ceiling(xmax/square_size)
      ystep = ceiling(ymax/square_size)
      kde_grid_size = c(ystep,xstep )
      
      
      W = owin(c(xmin,xmax),c(ymin,ymax))
      ref_ppp = as.ppp(ref,W)
      nonref_ppp = as.ppp(nonref,W)
      
      kde_ref = density(ref_ppp, positive = T, dimyx=kde_grid_size)
      plot(kde_ref,main="Ref Cell Intensity")
      contour(kde_ref)
      W_ref = levelset(kde_ref, thres, ">")
      
      kde_nonref = density(nonref_ppp, positive = T, dimyx=kde_grid_size)
      plot(kde_nonref,main="Nonref Cell Intensity")
      contour(kde_nonref,add=TRUE,col="blue")
      
      rownames(W_ref$m) = W_ref$yrow
      colnames(W_ref$m) = W_ref$xcol
      spdf_coord = coordinates(melt(W_ref$m))
      ind_inc = which(spdf_coord[,3]==1)
      spdf_coord = spdf_coord[ind_inc,-3]
      
      spdf_data = data.frame(ref = c(kde_ref$v), nonref = c(kde_nonref$v))
      spdf_data = spdf_data[ind_inc,]
      spdf_data$nonref[spdf_data$nonref<thres] = thres
      
      #TODO implement saving this out
      if(is.empty(spdf_coord[,1])){
        print("no common region")
        next()
      }
      
      spdf = SpatialPointsDataFrame(spdf_coord, spdf_data)
      
      dist.gw = gw.dist(spdf_coord)
      
      gwr.bw = bw.gwr(nonref~ref,
                      data = spdf,
                      approach = "AICc",
                      kernel = "gaussian",
                      adaptive = TRUE,
                      dMat = dist.gw)
      gwr_res = gwr.basic(nonref~0+ref, data = spdf, bw=gwr.bw, dMat = dist.gw)
      plot(density(gwr_res$SDF$ref))
      
      combined <- as.data.frame(cbind(spdf_coord, gwr_res$SDF$ref))
      colnames(combined) <- c("X_min", "Y_min", "Density")
      combined$X_max <- combined$X_min + xstep 
      combined$Y_max <- combined$Y_min + ystep
      ggplot() + geom_tile(aes(x=Y_min, y=X_min, fill=Density), data=combined) + geom_point(data=analysis_data, mapping=aes(x=Cell.X.Position, y=Cell.Y.Position))
      
      write.table(combined, paste(output_dir, paste(group_name, patient_id, nonrefcell, sep="_"), ".csv", sep=""), quote=FALSE, row.names=FALSE, sep=",")
      
    }
    
    
    
  }
}

# ggplot() + geom_point(aes(x=Var1, y=Var2, fill=V3), data=combined)
# ggplot() + geom_tile(aes(x=Var1, y=Var2, fill=V3), data=combined)

# DenMatrix =append(file,gwr_res$SDF$tum)
# IntMatrix=append(file,gwr_res$SDF$Intercept)
# a=list(tum = DenMatrix,
# int = IntMatrix)
# FinalOutput=c(FinalOutput,a)
