###############################################
# EXAMPLE:

setwd("/home/claudia/Documenti/lavori/Metilazione/ATLAS_epigenetic/Normalization/regRCPqn/R")
source("regRCPqn_regRCPqnREF.R")

# Example file containing M-values of 20 samples from GSE111629 pre-processed with rcp
fileName_M_data2 <- "../example/M_data2_part.txt"
M_data <- data.frame(fread(fileName_M_data2))
ref_path = "../example/"
data_name <- "Example"
# Run regRCPqnREF
M_data_norm_ref <- regRCPqnREF(M_data, ref_path, data_name)
###############################################
