# regRCPqn

An R package that extends RCP for meta-analysis usage introducing:

- a genomic region based normalization
- a between-samples normalization
- the possibility to normalize samples based on a reference distribution

## Installation:

Depends: data.table, preprocessCore, IlluminaHumanMethylation450kanno.ilmn12.hg19

library(devtools)
install_github("regRCPqn")

## regRCPqn:

regRCPqn(M_data, ref_path, data_name, save_ref=TRUE)

regRCPqnREF(M_data, ref_path, data_name)

### INPUT DESCRIPTION:

- M_data: data.frame containing M-values. The first column of M_data has to be named ID_REF and to contain the CpG codes. Each other column corresponds to a sample.

- ref_path: path to folder in which files will be saved

- data_name: prefix of the reference distribuion file name.

- save_ref (TRUE): whether to save the reference distribution of the normalized data set.


## Examples:

Usage example of regRCPqn and regRCPqnREF can be find in:

### EXAMPLE - regRCPqnREF

```
# Load data
data(BloodParkinson1)

# Set directory to save reference distribution
ref_path <- "./example/"
# Set dataset label
data_name <- "Example"

# Run regRCPqn
M_data_norm <- regRCPqn(M_data=BloodParkinson1, ref_path=ref_path, data_name=data_name, save_ref=TRUE)
```

### EXAMPLE - regRCPqnREF:

```
# Load data
data(BloodParkinson2)

# Set directory to save reference distribution
ref_path <- "./example/"
# Set dataset label
data_name <- "Example"

# Run regRCPqnREF
M_data_norm_ref <- regRCPqnREF(M_data=BloodParkinson2, ref_path=ref_path, data_name=data_name)
```

Examples include samples from Chuang YH, Paul KC, Bronstein JM, Bordelon Y et al. Parkinson's disease is associated with DNA methylation levels in human blood and saliva. Genome Med 2017 Aug 30;9(1):76. (GEO accession number GSE111629). Only a subset of CpG probes is included in the test files.

