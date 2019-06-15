###############################################
# EXAMPLE:

setwd("/home/claudia/Documenti/lavori/Metilazione/ATLAS_epigenetic/Normalization/regRCPqn/R")
source("regRCPqn_regRCPqnREF.R")

# Example file containing M-values of 20 samples from GSE111629 pre-processed with swan
fileName_M_data1 <- "../example/M_data1.txt"
M_data <- data.frame(fread(fileName_M_data1))
ref_path <- "../example/"
data_name <- "Example"
# Run regRCPqn
M_data_norm <- regRCPqn(M_data=M_data, ref_path=ref_path, data_name=data_name, save_ref=TRUE)
###############################################
