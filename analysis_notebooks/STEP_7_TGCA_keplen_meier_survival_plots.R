install.packages("survminer")
install.packages("survival")
BiocManager::install("DESeq2")



# script to run survival analysis using TCGA data
# setwd("~/Desktop/demo/survivalAnalysis")


library(TCGAbiolinks)
library(survminer)
library(survival)
library(SummarizedExperiment)
library(tidyverse)
library(DESeq2)

# getting clinical data for TCGA-BRCA cohort -------------------
clinical_brca <- GDCquery_clinic("TCGA-GBM")
any(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_brca) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_brca[,c(9,106,111)]


# looking at some variables associated with survival 

table(clinical_brca$vital_status)
# days_to_death, that is the number of days passed from the initial diagnosis to the patient's death (clearly, this is only relevant for dead patients)
# days_to_last_follow_up that is the number of days passed from the initial diagnosis to the last visit.

# change certain values the way they are encoded
clinical_brca$deceased <- ifelse(clinical_brca$vital_status == "Alive", FALSE, TRUE)

# create an "overall survival" variable that is equal to days_to_death
# for dead patients, and to days_to_last_follow_up for patients who
# are still alive
clinical_brca$overall_survival <- ifelse(clinical_brca$vital_status == "Alive",
                                         clinical_brca$days_to_last_follow_up,
                                         clinical_brca$days_to_death)





# get gene expression data -----------

# build a query to get gene expression data for entire cohort
query_brca_all = GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  #sample.type = "Primary Tumor",
  access = "open")

output_brca <- getResults(query_brca_all)
# get 20 primary tissue sample barcodes
tumor <- output_brca$cases#[1:20]
# OR
tumor <- output_brca[output_brca$sample_type == "Primary Tumor", "cases"]
tumor

# # get gene expression data from 20 primary tumors 
query <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling", # parameter enforced by GDCquery
  experimental.strategy = "RNA-Seq",
  workflow.type = "STAR - Counts",
  data.type = "Gene Expression Quantification",
  sample.type = c("Primary Tumor"),
  access = "open",
  barcode = tumor)

# download data
GDCdownload(query_brca)

########new code"

metadata <- query[[1]][[1]]

# Download data using api
GDCdownload(query, method = "api")

project_name='C:/Users/sergio.salas/Documents/GDCdata/TCGA-GBM/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification'

# Get main directory where data is stored
main_dir <- file.path(project_name)
# Get file list of downloaded files
file_list <- file.path(main_dir,list.files(main_dir,recursive = TRUE))
samples<-list.files(main_dir,recursive = FALSE)

# Read first downloaded to get gene names
test_tab <- read.table(file = file_list[1], sep = '\t', header = TRUE)
# Delete header lines that don't contain usefull information
test_tab <- test_tab[-c(1:4),]
# STAR counts and tpm datasets
tpm_data_frame <- data.frame(test_tab[,1])
count_data_frame <- data.frame(test_tab[,1])
tpkm_data_frame <- data.frame(test_tab[,1])

# Append cycle to get the complete matrix
for (i in c(1:length(file_list))) {
  # Read table
  test_tab <- read.table(file = file_list[i], sep = '\t', header = TRUE)
  # Delete not useful lines
  test_tab <- test_tab[-c(1:4),]
  # Column bind of tpm and counts data
  tpm_data_frame <- cbind(tpm_data_frame, test_tab[,7])
  tpkm_data_frame <- cbind(tpkm_data_frame, test_tab[,8])
  count_data_frame <- cbind(count_data_frame, test_tab[,4])
  # Print progres from 0 to 1
  print(i/length(file_list))
}





#############
# get counts
#tcga_brca_data <- GDCprepare(query_brca, summarizedExperiment = True)
brca_matrix <- tpkm_data_frame[,-1]#assay(tcga_brca_data, "unstranded")



# extract gene and sample metadata from summarizedExperiment object
gene_metadata<-as.data.frame(test_tab[,c('gene_id','gene_name')])
#gene_metadata <- as.data.frame(rowData(tcga_brca_data))
coldata <- as.data.frame(samples)
colnames(brca_matrix)=samples
rownames(brca_matrix)=test_tab[,'gene_id']

# vst transform counts to be used in survival analysis ---------------
# Setting up countData object   
dds <- DESeqDataSetFromMatrix(countData = brca_matrix,
                              colData = coldata,
                              design = ~ 1)

# Removing genes with sum total of 10 reads across all samples
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]


# vst 
vsd <- vst(dds, blind=FALSE)
brca_matrix_vst <- assay(vsd)
brca_matrix_vst[1:10,1:10]



# Get data for TP53 gene and add gene metadata information to it -------------
brca_tp53 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "NTRK2")

VEGFA <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "VEGFA")

HK2 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "HK2")

COL6A1 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "COL6A1")

MCAM <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "MCAM")


NOTCH3 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "NOTCH3")
GFAP <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "GFAP")
AQP4 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "AQP4")

TLR4 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "TLR4")


MS4A1 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "MS4A1")


SPARC <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "SPARC")

VIM <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "VIM")

MMP14 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "MMP14")


CXCL8 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "CXCL8")


LGALS3 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "LGALS3")


COL6A2 <- brca_matrix_vst %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "COL6A2")




brca_tp53$NTRK2<-brca_tp53$counts
brca_tp53$COL6A1<-COL6A1$counts
brca_tp53$VEGFA<-VEGFA$counts
brca_tp53$HK2<-HK2$counts
brca_tp53$MCAM<-MCAM$counts
brca_tp53$NOTCH3<-NOTCH3$counts
brca_tp53$GFAP<-GFAP$counts
brca_tp53$AQP4<-AQP4$counts
brca_tp53$MS4A1<-MS4A1$counts
brca_tp53$TLR4<-TLR4$counts
brca_tp53$SPARC<-SPARC$counts
brca_tp53$VIM<-VIM$counts
brca_tp53$MMP14<-MMP14$counts
brca_tp53$CXCL8<-CXCL8$counts
brca_tp53$LGALS3<-LGALS3$counts
brca_tp53$COL6A2<-COL6A2$counts



brca_tp53$module_expression <-rowMeans(brca_tp53[, c('GFAP','AQP4','NTRK2','MS4A1','TLR4')])

brca_tp53$module_expression <-rowMeans(brca_tp53[, c('COL6A1','SPARC','VIM')])

dat=brca_tp53[, c('COL6A1','VEGFA','HK2','MCAM','NOTCH3','MMP14','CXCL8','LGALS3')]
scaled.dat <- scale(dat)

# check that we get mean of 0 and sd of 1
colMeans(scaled.dat)  # faster version of apply(scaled.dat, 2, mean)
apply(scaled.dat, 2, sd)

brca_tp53[, c('COL6A1','VEGFA','HK2','MCAM','NOTCH3','MMP14','CXCL8','LGALS3')]<-scaled.dat

brca_tp53$module_expression <-rowMeans(brca_tp53[, c('COL6A1','VEGFA','HK2','MCAM','NOTCH3','MMP14','CXCL8','LGALS3')])



median_value <-quantile(brca_tp53$module_expression,probs=0.5)
#median_value_NTRK2 <-quantile(brca_tp53$NTRK2,probs=0.5)
perc75 <-quantile(brca_tp53$module_expression,probs=0.75)
perc25 <-quantile(brca_tp53$module_expression,probs=0.25)
brca_tp53$module_expression <- ifelse(brca_tp53$module_expression <= median_value,
                                      ifelse(brca_tp53$module_expression <= perc25,
                                             "LOW", "LOW"), ifelse(brca_tp53$module_expression <= perc75,
                                                     "HIGH", "HIGH"))

# Add clinical information to brca_tp53
brca_tp53$case_id <- gsub('-01.*', '', brca_tp53$case_id)
metadata_sample<-merge(brca_tp53, output_brca, by.x = 'case_id', by.y = 'id')
brca_tp53sub <- merge(metadata_sample, clinical_brca, by.x = 'cases.submitter_id', by.y = 'submitter_id')

dim(brca_tp53sub)


#write.csv(brca_tp53sub,'C:/Users/sergio.salas/Downloads/survival_data_gbm.csv')

# fitting survival curve -----------
fit <- survfit(Surv(overall_survival, deceased) ~ module_expression, data = brca_tp53sub)
fit
ggsurvplot(fit,
           data = brca_tp53sub,
           pval = T,
           risk.table = T)

### FITTING IT WITH MY DATA
##############3other stuff

# get median value
median_value <- median(brca_tp53$VEGFA)
brca_tp53$VEGFA <- ifelse(brca_tp53$VEGFA >= median_value, "HIGH", "LOW")
median_value <- median(brca_tp53$NTRK2)
brca_tp53$NTRK2 <- ifelse(brca_tp53$NTRK2 >= median_value, "HIGH", "LOW")
median_value <- median(brca_tp53$COL6A1)
brca_tp53$COL6A1 <- ifelse(brca_tp53$COL6A1 >= median_value, "HIGH", "LOW")







genes_of_interest <- c("", "NTRK2", "gene3")  # Add your gene names here

# Fit a multivariate Kaplan-Meier survival model
surv_fit <- survfit(Surv(overall_survival, deceased) ~ VEGFA + NTRK2+COL6A1, data = brca_tp53sub)
# Perform multivariate log-rank test
p_value_multivar <- survdiff(Surv(overall_survival, deceased) ~ VEGFA + NTRK2+COL6A1, data = brca_tp53sub)$chisq[1]

print(p_value_multivar)

ggsurvplot(
  surv_fit,
  data = brca_tp53sub,
  title = "Multivariate Survival Curves",
  xlab = "Time",
  ylab = "Survival Probability"
)




fit2 <- survdiff(Surv(overall_survival, deceased) ~ strata, data = brca_tp53)