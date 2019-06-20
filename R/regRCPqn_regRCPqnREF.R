###############################################
### INPUT DESCRIPTION:

# M_data: data.frame containing M-values. The first column of M_data has to be named ID_REF and to contain the CpG codes. Each other column corresponds to a sample.
# ref_path: path to folder in which files will be saved
# data_name: prefix that will be used to create the file name of the reference distribuion if save_ref is TRUE.
# save_ref (TRUE): whether to save the reference distribution of the normalized data set.

###############################################
# %%
library(data.table)
regRCPqn <- function(M_data, ref_path, data_name,save_ref=TRUE, save_norm=TRUE){
  ###############################################
  # Set default settings of RCP
  ###############################################
  dist_ <- 25
  quantile.grid <- seq(0.001,0.999,by=0.001)
  ###############################################
  # Read annotation file
  # Remove chromosomes X and Y
  # Select the same CpG in the annotation file and in the file with M-values
  ###############################################
  i <-1
  annot450k <- fread(paste0("../input/HumanMethylation450_15017482_v1-2_CLEAN_part",i,".csv"),sep="\t")
  for (i in 2:9){
    annot_part <- fread(paste0("../input/HumanMethylation450_15017482_v1-2_CLEAN_part",i,".csv"),sep="\t")
    annot450k <- rbind(annot450k,annot_part)
  }
  annot450k <- data.frame(annot450k)
  annot450k <- annot450k[! annot450k$CHR %in% c("X","Y"), ]
  rownames(annot450k) <- annot450k$IlmnID
  rownames(M_data) <- M_data$ID_REF
  cpg_kept <- intersect(rownames(annot450k), rownames(M_data))
  annot450k <- annot450k[cpg_kept,]
  annot450k <- annot450k[order(annot450k$CHR, annot450k$MAPINFO),]
  M_data <- M_data[rownames(annot450k),]
  M_data <- M_data[,colnames(M_data)!="ID_REF"]
  ###############################################
  # Select genomic region types
  ###############################################
  regionTypes <- as.vector(annot450k$Relation_to_UCSC_CpG_Island)
  regionTypes[regionTypes %in% c("N_Shore", "S_Shore")] <- "Shore"
  regionTypes[regionTypes %in% c("N_Shelf", "S_Shelf")] <- "Shelf"
  regionTypes[regionTypes==""] <- "Ocean"
  annot450k$Relation_to_UCSC_CpG_Island_Summary <- regionTypes
  regionTypes <- unique(regionTypes)
  annot450k$Relation_to_UCSC_CpG_Island_Summary
  ###############################################
  # Separately for each region type and quantile normalize samples:
  # - Run RCP:
  #     - Select a subset of close probes with different probe type.
  #     - Use the M-values of the selected probes to indentify the transformation between type II and type I probes.
  #     - Transform all type II probes according to that transformation.
  # - Separately for each probe type quantile normalize samples.
  ###############################################
  for (region_i in seq(length(regionTypes))){
    annot450k_part <- annot450k[annot450k$Relation_to_UCSC_CpG_Island_Summary == regionTypes[region_i], ]
    probe.II.Name <- annot450k_part$Name[annot450k_part$Infinium_Design_Type == "II"]
    probe.I.Name <- annot450k_part$Name[annot450k_part$Infinium_Design_Type == "I"]
    anno1 <- annot450k_part[1:(nrow(annot450k_part)-1),]
    anno2 <- annot450k_part[2:nrow(annot450k_part),]
    flag <- ((abs(anno1$MAPINFO - anno2$MAPINFO)) < dist_ & (anno1$CHR == anno2$CHR) &
     (anno1$Infinium_Design_Type != anno2$Infinium_Design_Type))
    anno1 <- anno1[flag,]
    anno2 <- anno2[flag,]
    probe.I <- rownames(anno1)
    probe.II <- rownames(anno2)
    probe.I[anno2$Infinium_Design_Type=="I"] <- rownames(anno2)[anno2$Infinium_Design_Type == "I"]
    probe.II[anno1$Infinium_Design_Type=="II"] <- rownames(anno1)[anno1$Infinium_Design_Type == "II"]
    raw.M.t <- M_data[c(probe.I, probe.II),]
    #linear regression
    M.II <- raw.M.t[probe.II,]
    M.I <- raw.M.t[probe.I,]
    qtl <- function(x) quantile(x, quantile.grid, na.rm=TRUE)
    M.I <- apply(M.I,2,qtl)  # compute quantiles for M.I
    M.II <- apply(M.II,2,qtl)  # compute quantiles for M.II
    beta.est <- mat.or.vec(2,ncol(M_data))
    # for each sample select the finite quantiles of type I (Y) and of type II (X)
    # and estimate the regression coefficient of the model Y = beta.est * X
    for (i in 1:ncol(M_data)){
      index <- ((M.II[,i]!=Inf) & (M.II[,i]!=-Inf) & (M.I[,i]!=Inf) & (M.I[,i]!=-Inf))
      X <- cbind(rep(1, sum(index)), M.II[index, i])
      Y <- M.I[index, i]
      beta.est[, i] <- solve(t(X)%*%X)%*%t(X)%*%Y
    }
    # Separately for each sample, map all type II probes with the fitted transformation
    M.II.all <- M_data[which(rownames(M_data) %in% probe.II.Name),]
    M.II.new <- mat.or.vec(nrow(M.II.all),ncol(M.II.all))
    for (i in 1:ncol(M.II.all)){
      M.II.new[,i] <- beta.est[1,i] + beta.est[2,i] * M.II.all[,i]
    }
    M.II.new[M.II.all == Inf] <- Inf
    M.II.new[M.II.all == -Inf] <- (-Inf)
    # Replace M-values of the type II probes of the current genomic region type with the transformed data.
    M_data[which(rownames(M_data) %in% probe.II.Name),] <- M.II.new
    # Quantile normalization between samples computed separately for type I and type II probes
    M_data[which(rownames(M_data) %in% probe.II.Name),] <- preprocessCore::normalize.quantiles(as.matrix(M_data[which(rownames(M_data) %in% probe.II.Name),]))
    M_data[which(rownames(M_data) %in% probe.I.Name),] <- preprocessCore::normalize.quantiles(as.matrix(M_data[which(rownames(M_data) %in% probe.I.Name),]))
    # Compute and save target distributions:
    if(save_ref){
      target_distr_II <- rowMeans(M_data[which(rownames(M_data) %in% probe.II.Name),])
      target_distr_II <- data.frame(ID_REF=names(target_distr_II), values=as.vector(target_distr_II))
      target_distr_I <- rowMeans(M_data[which(rownames(M_data) %in% probe.I.Name),])
      target_distr_I <- data.frame(ID_REF=names(target_distr_I), values=as.vector(target_distr_I))
      fwrite(target_distr_II, file=paste(ref_path,data_name,"_target_distribution_typeII_",regionTypes[region_i],".txt",sep=""))
      fwrite(target_distr_I, file=paste(ref_path,data_name,"_target_distribution_typeI_",regionTypes[region_i],".txt",sep=""))
    }
  }
  return(M_data)

}

# %%
###############################################
# INPUT DESCRIPTION:

# M_data: data.frame containing M-values. The first column of M_data has to be named ID_REF and to contain the CpG codes. Each other column corresponds to a sample.
# ref_path: path to folder in which the reference distribution of the genomic region types computed with regRCPqn have been saved.
# data_name: prefix used in regRCPqn to save the reference distribution files

###############################################
regRCPqnREF <- function(M_data, ref_path, data_name){
  ###############################################
  # Set default settings of RCP
  ###############################################
  dist_ <- 25
  quantile.grid <- seq(0.001,0.999,by=0.001)
  ###############################################
  # Read annotation file
  # Remove chromosomes X and Y
  # Select the same CpG in the annotation file and in the file with M-values
  ###############################################
  i <-1
  annot450k <- fread(paste0("../input/HumanMethylation450_15017482_v1-2_CLEAN_part",i,".csv"),sep="\t")
  for (i in 2:9){
    annot_part <- fread(paste0("../input/HumanMethylation450_15017482_v1-2_CLEAN_part",i,".csv"),sep="\t")
    annot450k <- rbind(annot450k,annot_part)
  }
  annot450k <- data.frame(annot450k)
  annot450k <- annot450k[! annot450k$CHR %in% c("X","Y"), ]
  rownames(annot450k) <- annot450k$IlmnID
  rownames(M_data) <- M_data$ID_REF
  cpg_kept <- intersect(rownames(annot450k), rownames(M_data))
  annot450k <- annot450k[cpg_kept,]
  annot450k <- annot450k[order(annot450k$CHR, annot450k$MAPINFO),]
  M_data <- M_data[rownames(annot450k),]
  M_data <- M_data[,colnames(M_data)!="ID_REF"]
  ###############################################
  # Select genomic region types
  ###############################################
  regionTypes <- as.vector(annot450k$Relation_to_UCSC_CpG_Island)
  regionTypes[regionTypes %in% c("N_Shore", "S_Shore")] <- "Shore"
  regionTypes[regionTypes %in% c("N_Shelf", "S_Shelf")] <- "Shelf"
  regionTypes[regionTypes==""] <- "Ocean"
  annot450k$Relation_to_UCSC_CpG_Island_Summary <- regionTypes
  regionTypes <- unique(regionTypes)
  ###############################################
  # Separately for each region type and quantile normalize samples:
  # - Run RCP:
  #     - Select a subset of close probes with different probe type.
  #     - Use the M-values of the selected probes to indentify the transformation between type II and type I probes.
  #     - Transform all type II probes according to that transformation.
  # - Separately for each probe type quantile normalize samples on target distribution.
  ###############################################
  for (region_i in seq(length(regionTypes))){
    annot450k_part <- annot450k[annot450k$Relation_to_UCSC_CpG_Island_Summary == regionTypes[region_i], ]
    probe.II.Name <- annot450k_part$Name[annot450k_part$Infinium_Design_Type == "II"]
    probe.I.Name <- annot450k_part$Name[annot450k_part$Infinium_Design_Type == "I"]
    anno1 <- annot450k_part[1:(nrow(annot450k_part)-1),]
    anno2 <- annot450k_part[2:nrow(annot450k_part),]
    flag <- ((abs(anno1$MAPINFO - anno2$MAPINFO)) < dist_ & (anno1$CHR == anno2$CHR) &
     (anno1$Infinium_Design_Type != anno2$Infinium_Design_Type))
    anno1 <- anno1[flag,]
    anno2 <- anno2[flag,]
    probe.I <- rownames(anno1)
    probe.II <- rownames(anno2)
    probe.I[anno2$Infinium_Design_Type=="I"] <- rownames(anno2)[anno2$Infinium_Design_Type == "I"]
    probe.II[anno1$Infinium_Design_Type=="II"] <- rownames(anno1)[anno1$Infinium_Design_Type == "II"]
    raw.M.t <- M_data[c(probe.I, probe.II),]
    # linear regression
    M.II <- raw.M.t[probe.II,]
    M.I <- raw.M.t[probe.I,]
    qtl <- function(x) quantile(x, quantile.grid, na.rm=TRUE)
    M.I <- apply(M.I,2,qtl)  # compute quantiles for M.I
    M.II <- apply(M.II,2,qtl)    # compute quantiles for M.II
    beta.est <- mat.or.vec(2,ncol(M_data))
    # for each sample select the finite quantiles of type I (Y) and of type II (X)
    # and estimate the regression coefficient of the model Y = beta.est * X
    for (i in 1:ncol(M_data)){
      index <- ((M.II[,i]!=Inf) & (M.II[,i]!=-Inf) & (M.I[,i]!=Inf) & (M.I[,i]!=-Inf))
      X <- cbind(rep(1, sum(index)), M.II[index, i])
      Y <- M.I[index, i]
      beta.est[, i] <- solve(t(X)%*%X)%*%t(X)%*%Y
    }
    # Separately for each sample, map all type II probes with the fitted transformation
    M.II.all <- M_data[which(rownames(M_data) %in% probe.II.Name),]
    M.II.new <- mat.or.vec(nrow(M.II.all),ncol(M.II.all))
    for (i in 1:ncol(M.II.all)){
      M.II.new[,i] <- beta.est[1,i] + beta.est[2,i] * M.II.all[,i]
    }
    M.II.new[M.II.all == Inf] <- Inf
    M.II.new[M.II.all == -Inf] <- (-Inf)
    # Replace M-values of the type II probes of the current genomic region type with the transformed data
    M_data[which(rownames(M_data) %in% probe.II.Name),] <- M.II.new
    # Quantile normalize samples on reference distribution separately for probe type
    ref_typeII <- paste(ref_path,data_name,"_target_distribution_typeII_",regionTypes[region_i],".txt",sep="")
    ref_typeI <- paste(ref_path,data_name,"_target_distribution_typeI_",regionTypes[region_i],".txt",sep="")
    target_distr_II <- data.frame(fread(ref_typeII))
    target_distr_I <- data.frame(fread(ref_typeI))
    target_distr_II <- target_distr_II[which(target_distr_II$ID_REF %in% rownames(annot450k)),]
    target_distr_I <- target_distr_I[which(target_distr_I$ID_REF %in% rownames(annot450k)),]
    M_data[which(rownames(M_data) %in% probe.II.Name),] <- preprocessCore::normalize.quantiles.use.target(as.matrix(M_data[which(rownames(M_data) %in% probe.II.Name),]),target = target_distr_II$value)
    M_data[which(rownames(M_data) %in% probe.I.Name),] <- preprocessCore::normalize.quantiles.use.target(as.matrix(M_data[which(rownames(M_data) %in% probe.I.Name),]),target = target_distr_I$value)
  }
  return(M_data)
}

# %%
