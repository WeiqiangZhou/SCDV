# SCDV: Single cell differential variance analysis

### Overview
SCDV is an R package for analyzing differentially variable genes in single-cell RNA-seq data.

### Installation
To install SCDV via Github, run the following code in R: 
```
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("WeiqiangZhou/SCDV")
```

### Updates
Updated **differential mean test** and **differential var test** by using the log-transformed data in the tests.

### How to use
There are three major steps for using SCDV

#### Step 1. Estimate dropout probability
```
library(SCDV)

##input your count data as data_treatment_count and data_control_count, rownames of the input matrix should be ensembl id
data_treatment_count <- input_data_treatment
data_control_count <- input_data_control

##get gene information
gene_info <- annotation_data[["gene_len_GRCh37"]] 
##change to "gene_info <- annotation_data[["gene_len_GRCm38"]]" for mouse data

match_gene_idx <- match(rownames(data_treatment_count),gene_info[,1])
gene_len <- gene_info[match_gene_idx,2]
names(gene_len) <- gene_info[match_gene_idx,1]
match_gene_name <- gene_info[match_gene_idx,-2]

##get neighboring cells
neighbor_result_treatment <- get_expected_cell(data_treatment_count,gene_len)
neighbor_result_control <- get_expected_cell(data_control_count,gene_len)

##estimate dropout probability, set multi-core using ncore, if not using multi-core, set ncore=1
output_treatment <- estimate_dropout_main(data_treatment_count,neighbor_result_treatment$data_expect,gene_len,ncore=6)
output_control <- estimate_dropout_main(data_control_count,neighbor_result_control$data_expect,gene_len,ncore=6)

save(output_treatment,file="output_treatment.rda")
save(output_control,file="output_control.rda")
```

#### Step 2. Adjust library size using house keeping genes

```
hp_gene <- annotation_data[["LV_gene_human"]][,1]
##change to "hp_gene <- annotation_data[["LV_gene_mouse"]][,1]" for mouse data

gene_names <- sapply(rownames(data_treatment_count),function(x) sub("\\..*","",x))
library_scale <- adjust_library_size(c(output_treatment,output_control),hp_gene,gene_names)

treatment_data_true <- sapply(output_treatment,function(x) x$data_true)
treatment_data_weight <- 1 - sapply(output_treatment,function(x) x$post_weight)

control_data_true <- sapply(output_control,function(x) x$data_true)
control_data_weight <- 1 - sapply(output_control,function(x) x$post_weight)

treatment_data_adjust <- t(t(treatment_data_true)*library_scale[1:length(output_treatment)])
control_data_adjust <- t(t(control_data_true)*library_scale[-c(1:length(output_treatment))])
```

#### Step 3. Perform differential analysis
SCDV supports three types of tests: 

#### 1. differential hyper-variability test
```
##set multi-core using ncore, if not using multi-core, set ncore=1
diff_disper <- scdv_main(treatment_data_adjust,treatment_data_weight,control_data_adjust,control_data_weight,num_permute=10000,span_param=0.5,ncore=6)
write.csv(cbind(match_gene_name,diff_disper),file="diff_hypervar.csv",row.names=FALSE)
```

#### 2. differential variability test
```
##set multi-core using ncore, if not using multi-core, set ncore=1
diff_var <- test_var_main(treatment_data_adjust,treatment_data_weight,control_data_adjust,control_data_weight,num_permute=10000,ncore=6,log_transform=TRUE)
write.csv(data.frame(match_gene_name,diff_var),file="diff_var.csv",row.names=FALSE)
```

#### 3. differential mean test
```
##set multi-core using ncore, if not using multi-core, set ncore=1
diff_expr <- test_mean_main(treatment_data_adjust,treatment_data_weight,control_data_adjust,control_data_weight,num_permute=10000,ncore=6,log_transform=TRUE)
write.csv(data.frame(match_gene_name,diff_expr),file="diff_expr.csv",row.names=FALSE)
```

### Q&A
#### What's the assumption of the tests in SCDV?
All the tests in SCDV are permutation-based tests. These are non-parametric tests that do not rely on strong parametric assumptions.

#### What's the difference between hyper-variability test and variability test?
The major difference is they test different types of variance.

The variance of a gene in single-cell RNA-seq data has been showed to be highly correlated with the mean expression of the gene. To remove such "mean effects", SCDV calculates the hyper-variability statistics which is a ratio of the observed variance of a gene to the expected variance at the same mean expression level.

#### What are the return values of each test?
Use `?scdv_main`, `?test_var_main`, and `?test_mean_main` in R to check the description of the return values in the help page of the package.

#### What is "house keeping genes"?
The "house keeping genes" are genes that show low variability but consistently expressed in different cell types which are obtained using recount2 (https://jhubiostatistics.shinyapps.io/recount/).

#### How to cite SCDV?
Use `citation("SCDV")` in R to get the citation information.

### Contact
Please contact **Weiqiang Zhou**: wzhou14@jhu.edu for questions and suggestions.
