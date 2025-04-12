#### pkgs ####
require(data.table)
require(glue)
require(DESeq2)
require(readxl)
require(ggvenn)
require(readxl)
require(limma)
require(oligo)
require(pd.hugene.1.0.st.v1)
require(biomaRt)
require(GEOquery)
require(tidyverse)

expr_df <- list()
meta_df <- list()
panel_df <- list()
str_df <- list()

# GSE32472
## chip; log transformed
str_df[["gse32472"]] <- "chip; log transformed"

panel_df[["gse32472"]] <- read.delim("GEO_data/GSE32472/GPL6244-17930.txt",
                                     header = TRUE, sep = "\t",
                                     comment.char = "#"
) %>%
  filter(gene_assignment != "---") %>%
  separate(gene_assignment,
           c("GB", "gene_name", "category"),
           sep = " // "
  ) %>%
  select(ID, gene_name)

meta_df[["gse32472"]] <- fread("GEO_data/GSE32472/meta_data_gse32472.csv") %>%
  na.omit() %>%
  separate(`bronchopulmonary dysplasia (bpd) group:ch1`,
           into = c("BPD_status", "BPD_type"),
           sep = "\\s*\\(|\\)"
  ) %>%
  rename(gse_id = "V1") %>%
  select(gse_id, BPD_type) %>%
  na.omit() %>%
  mutate(
    BPD_status = ifelse(BPD_type == "no BPD", "non_BPD", "BPD"),
    BPD_type = str_replace(BPD_type, " ", "_"),
    BPD_type = str_replace(BPD_type, "no_", "non_")
  )

cel_path <- "GEO_data/GSE32472/"
cel_files <- list.files(cel_path, pattern = "\\.CEL.gz$", full.names = TRUE)
GSE32472 <- read.celfiles(cel_files, pkgname = "pd.hugene.1.0.st.v1")
GSE32472_rma <- rma(GSE32472)
expr_cel_rma <- exprs(GSE32472_rma)
colnames(expr_cel_rma) <- gsub("^(GSM\\d+)_.*$", "\\1", colnames(expr_cel_rma))

expr_df[["gse32472"]] <- as.data.frame(expr_cel_rma) %>%
  rownames_to_column("ID_REF") %>%
  mutate(ID_REF = as.numeric(ID_REF)) %>%
  right_join(panel_df[["gse32472"]], by = c("ID_REF" = "ID")) %>%
  select(-ID_REF) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  column_to_rownames("gene_name") %>%
  select(all_of(meta_df[["gse32472"]]$gse_id))

# GSE220135
## RNA-seq; raw counts
str_df[["gse220135"]] <- "RNA-seq; raw counts"

panel_df[["gse220135"]] <- fread("GEO_data/GSE220135/Human.GRCh38.p13.annot.tsv.gz") %>%
  select(GeneID, Symbol) %>%
  rename(
    ID = "GeneID",
    gene_name = Symbol
  )

meta_df[["gse220135"]] <- read_xlsx("GEO_data/GSE220135/GSE220135_sample information.xlsx") %>%
  select(Accession, `Variable-disease`, `Variable-time`) %>%
  rename(gse_id = "Accession", BPD_status = `Variable-disease`, time = `Variable-time`) %>%
  mutate(
    BPD_status = ifelse(BPD_status == "BPD", "BPD", "non_BPD"),
    time = ifelse(time == "CordBlood", "Day0", time)
  )

expr_df[["gse220135"]] <- fread("GEO_data/GSE220135/GSE220135_raw_counts_GRCh38.p13_NCBI.tsv") %>%
  left_join(panel_df[["gse220135"]], by = c("GeneID" = "ID")) %>%
  select(-GeneID) %>%
  group_by(gene_name) %>%
  summarise(across(everything(), mean)) %>%
  ungroup() %>%
  column_to_rownames("gene_name") %>%
  select(any_of(meta_df[["gse220135"]]$gse_id))

write_rds(expr_df, "bkup_data/expression_matrix_original.rds")
write_rds(meta_df, "bkup_data/meta_data.rds")

