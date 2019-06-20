# regRCPqn

An R package that extends RCP for meta-analysis usage introducing:

- a genomic region based normalization
- a between-samples normalization
- the possibility to normalize samples based on a reference distribution

## Installation:

Depends: data.table

library(devtools)
install_github("regRCPqn")

## regRCPqn:

regRCPqn(M_data, fileName_annot450k, ref_path, data_name,save_ref=TRUE, save_norm=TRUE)

regRCPqnREF(M_data, fileName_annot450k, ref_path, data_name)

### INPUT DESCRIPTION:

- M_data: data.frame containing M-values. The first column of M_data has to be named ID_REF and to contain the CpG codes. Each other column corresponds to a sample.

- ref_path: path to folder in which files will be saved

- data_name: prefix of the reference distribuion file name.

- save_ref (TRUE): whether to save the reference distribution of the normalized data set.


## Examples:

Usage example of regRCPqn and regRCPqnREF can be find in:

./R/test_regRCPqn.R

./R/test_regRCPqnREF.R

Examples include samples from Chuang YH, Paul KC, Bronstein JM, Bordelon Y et al. Parkinson's disease is associated with DNA methylation levels in human blood and saliva. Genome Med 2017 Aug 30;9(1):76. (GEO accession number GSE111629). Only a subset of CpG probes is included in the test files.

