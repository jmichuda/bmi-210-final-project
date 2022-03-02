library(maftools)
library(tidyverse)


# INPUT/OUTPUT FILES -----------------------------------------------------------
cosmic_cancer_gene_tbl <- "src/COSMIC_Census_allTue_Feb_22_00_32_56_2022.tsv"
oncokb_cancer_gene_tbl <- "src/OncoKB_cancerGeneList.tsv"
tcga_dna_repair_gene_tbl <- "src/TCGA_DNA_Repair_Genes_PMID_29617664.xlsx"
tcga_gcop_file = snakemake@input[["copyalt"]]
tcga_gene_fusion_tbl <- "src/TCGA_DriverFusions_mmc2.xlsx"

output_file = snakemake@output[["tcga_var_tbl"]]
output_file_clinical = snakemake@output[["tcga_clinical"]]
output_file_cna = snakemake@output[["tcga_cna"]]
output_file_fusions = snakemake@output[["tcga_fusions"]]
output_file_maf = snakemake@output[["tcga_maf"]]

# FILTERING CRITERIA -----------------------------------------------------------
RESTRICT_TO_CANCER_GENES = TRUE
# Somatic SNV/Indel
MIN_NORMAL_WXS_DEPTH <- 20
MIN_TUMOR_WXS_DEPTH <- 50
MIN_TUMOR_ALT_READ_COUNT <- 5
# Somatic Copy Number
FILTER_LOW_LEVEL_COPY_GAINS <- FALSE
# Gene Fusion
MIN_SPANNING_READ_COUNT <- 3
MIN_JUNCTION_READ_COUNT <- 2

# ASSEMBLE CANCER GENE LIST ----------------------------------------------------
message("---- Assembling Cancer Gene List ...")
cosmic_cangenes <- readr::read_tsv(cosmic_cancer_gene_tbl, show_col_types = FALSE)
oncokb_cangenes <- readr::read_tsv(oncokb_cancer_gene_tbl, show_col_types = FALSE)
tcga_dna_repair_genes <- readxl::read_xlsx(tcga_dna_repair_gene_tbl, skip = 3)

cangenes <- c(
  cosmic_cangenes$`Gene Symbol`,
  oncokb_cangenes$`Hugo Symbol`,
  tcga_dna_repair_genes$`Gene Symbol`
) %>%
  unique() %>%
  stringr::str_trim() %>%
  sort()

# UTILS ------------------------------------------------------------------------
# https://docs.gdc.cancer.gov/Encyclopedia/pages/TCGA_Barcode/
# Tumor_Sample_Barcode: TCGA-OR-A5J1-01A-11D-A29I-10
# Participant_ID: TCGA-OR-A5J1
# Sample_ID: TCGA-OR-A5J1-01

barcode_to_participant <- function(x) {
  stringr::str_split(x, pattern = "-") %>%
    purrr::map_chr(function(parts) {
      paste0(parts[1:3], collapse = "-")
    })
}


barcode_to_sample <- function(x) {
  stringr::str_split(x, pattern = "-") %>%
    purrr::map_chr(function(parts) {
      paste0(parts[1:4], collapse = "-")
    }) %>%
    purrr::map_chr(function(s) {
      stringr::str_sub(s, 1, nchar(s) - 1)
    })
}

# CONSTRUCT: TCGA PANCAN MAF (SNV/INDEL) DATA TABLE ----------------------------
message("---- Loading TCGA PanCan MAF (SNV/Indel) Data ...")
tcga_cohorts = maftools::tcgaAvailable() %>%
  filter(Study_Abbreviation != "Unknown") 

# Construct: tibble data frame with each maf column entry containing TCGA study 
# (cancer type) specific (non-silent) mutation data in maf format.
tcga_maf = tcga_cohorts %>%
  select(Study_Abbreviation, Study_Name) %>%
  mutate(
    maf = purrr::map(Study_Abbreviation, function(cohort) {
      maftools::tcgaLoad(cohort, source = "MC3")@data %>% 
        as_tibble() %>%
        mutate_if(is.factor, as.character)
    })
  )

# Assert: each maf table has the same (ordered) columns.
expected_cols <- colnames(tcga_maf$maf[[1]])
stopifnot(
  all(
    purrr::map_lgl(tcga_maf$maf, function(df) {
      identical(colnames(df), expected_cols)
    })
  )
)

tcga_maf <- tcga_maf %>% tidyr::unnest(cols = "maf")
#stopifnot(identical(tcga_maf$Study_Abbreviation, tcga_maf$Cohort))
#tcga_maf <- tcga_maf %>% select(-Cohort)

# FILTER SNV/INDELS ------------------------------------------------------------
message("---- Preparing SNV/Indel Data ...")
if (RESTRICT_TO_CANCER_GENES) {
  tcga_maf <- tcga_maf %>% filter(Hugo_Symbol %in% cangenes)
}

tcga_maf <- tcga_maf %>%
  filter(n_depth >= MIN_NORMAL_WXS_DEPTH) %>%
  filter(t_depth >= MIN_TUMOR_WXS_DEPTH)  %>%
  filter(t_alt_count >= MIN_TUMOR_ALT_READ_COUNT)

# -------------------------------------------------------------------
tcga_maf_out <- tcga_maf %>%
  mutate(
    # Set to sample-level id shared by copy and fusion data samples
    Tumor_Sample_Barcode = barcode_to_sample(Tumor_Sample_Barcode)
  )
readr::write_tsv(tcga_maf_out, file = gzfile(output_file_maf))

tcga_clinical <- tcga_maf_out %>%
  select(Tumor_Sample_Barcode, Study_Abbreviation) %>%
  rename(
    SAMPLE_ID = Tumor_Sample_Barcode,
    ONCOTREE_CODE = Study_Abbreviation
  ) %>%
  unique() %>%
  arrange(ONCOTREE_CODE, SAMPLE_ID)
readr::write_tsv(tcga_clinical, file = output_file_clinical)
# -------------------------------------------------------------------

tcga_maf <- tcga_maf %>%
  mutate(
    participant_id = barcode_to_participant(Tumor_Sample_Barcode),
    sample_id = barcode_to_sample(Tumor_Sample_Barcode)
  ) %>%
  select(participant_id, sample_id, Tumor_Sample_Barcode, everything())


# CONSTRUCT: TCGA PANCAN GENE-LEVEL COPY STATE TABLE ---------------------------
# Source:
# https://xenabrowser.net/datapages/?dataset=TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_thresholded.by_genes&host=https%3A%2F%2Ftcga.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

message("---- Preparing Gene-Level Copy Alteration Data ...")
tcga_mut_samples <- tcga_maf$sample_id %>% unique()

Sys.setenv(VROOM_CONNECTION_SIZE = "10000000 KiB")
tcga_cop <- readr::read_tsv(tcga_gcop_file, show_col_types = FALSE) %>%
  rename(gene = Sample) %>%
  filter(gene %in% cangenes)
# -------------------------------------------------------------------
tcga_cna_out <- tcga_cop %>%
  select(c("gene", intersect(colnames(.), tcga_mut_samples))) %>%
  rename(`Gene Symbol` = gene) %>%
  mutate(
    `Locus ID` = NA,
    Cytoband = NA
  ) %>%
  select(`Gene Symbol`, `Locus ID`, Cytoband, everything()) %>%
  arrange(`Gene Symbol`)
readr::write_tsv(tcga_cna_out, file = gzfile(output_file_cna))
# -------------------------------------------------------------------

tcga_cop <- tcga_cop %>%
  tidyr::pivot_longer(cols = setdiff(colnames(.), "gene"), 
                      names_to = "sample_id", values_to = "copy_state") %>%
  filter(copy_state != 0) %>%
  filter(sample_id %in% tcga_mut_samples) %>%
  select(sample_id, gene, copy_state) %>%
  arrange(sample_id, gene)

if (FILTER_LOW_LEVEL_COPY_GAINS) {
  tcga_cop <- filter(tcga_cop, tcga_cop$copy_state != 1)
}

# CONSTRUCT: TCGA PANCAN GENE FUSION DATA TABLE --------------------------------
# Source: Driver Fusions and Their Implications in the Development 
#           and Treatment of Human Cancers
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5916809/ (Table S1)

message("---- Preparing Gene Fusion Data ...")
tcga_fusion <- readxl::read_xlsx(tcga_gene_fusion_tbl, sheet = 2, skip = 1) %>%
  rename(
    Study_Abbreviation = Cancer,
    Tumor_Sample_Barcode = Sample
  ) %>%
  mutate(
    participant_id = barcode_to_participant(Tumor_Sample_Barcode),
    sample_id = barcode_to_sample(Tumor_Sample_Barcode),
    gene1 = purrr::map_chr(Fusion, function(s) {
      stringr::str_split(s, pattern = "--")[[1]][[1]]
    }),
    gene2 = purrr::map_chr(Fusion, function(s) {
      stringr::str_split(s, pattern = "--")[[1]][[2]]
    })
  ) %>%
  filter((gene1 %in% cangenes) | (gene2 %in% cangenes)) %>%
  filter(sample_id %in% tcga_mut_samples) %>%
  filter(Spanning >= MIN_SPANNING_READ_COUNT) %>%
  filter(Junction >= MIN_JUNCTION_READ_COUNT) %>%
  select(Study_Abbreviation, participant_id, sample_id, Tumor_Sample_Barcode,
         everything()) %>%
  arrange(Study_Abbreviation, participant_id, sample_id, Tumor_Sample_Barcode)

# -------------------------------------------------------------------
tcga_fusion_out <- tcga_fusion %>%
  select(sample_id, Fusion) %>%
  rename(Tumor_Sample_Barcode = sample_id) %>%
  mutate(
    Fusion = stringr::str_replace(Fusion, "--", "-")
  ) %>%
  arrange(Tumor_Sample_Barcode, Fusion)
readr::write_tsv(tcga_fusion_out, file = output_file_fusions)
# -------------------------------------------------------------------

# CONSTRUCT: TCGA PANCANCER COMBINED VARIANT DATA TABLE ------------------------
message("---- Constructing Combined Variant Data Table ...")
genvar_mut <- tcga_maf %>%
  select(Study_Abbreviation, participant_id, sample_id, Tumor_Sample_Barcode,
         Hugo_Symbol, Variant_Classification, HGVSp_Short) %>%
  rename(
    TCGA_Cohort = Study_Abbreviation, 
    Participant_ID = participant_id, 
    Sample_ID = sample_id, 
    Gene = Hugo_Symbol,
    Variant_Type = Variant_Classification,
    Description = HGVSp_Short
  )

sample_to_cohort <- genvar_mut %>%
  select(Participant_ID, Sample_ID, TCGA_Cohort) %>%
  unique() %>%
  arrange(TCGA_Cohort, Participant_ID, Sample_ID)

copy_state_to_description <- c(
  "-2" = "homozygous_deletion",
  "-1" = "single_copy_loss",
  "1"  = "low_copy_gain",
  "2"  = "high_copy_gain"
)

genvar_cop <- tcga_cop %>%
  left_join(sample_to_cohort, by = c("sample_id" = "Sample_ID")) %>%
  mutate(
    Tumor_Sample_Barcode = NA,
    Variant_Type = "Copy_Number_Alteration",
    Description = copy_state_to_description[as.character(copy_state)]
  ) %>%
  rename(
    Gene = gene,
    Sample_ID = sample_id
  ) %>%
  select(colnames(genvar_mut))

genvar_fus <- tcga_fusion %>%
  rename(
    TCGA_Cohort = Study_Abbreviation,
    Participant_ID = participant_id,
    Sample_ID = sample_id,
    Gene = Fusion
  ) %>%
  mutate(
    Variant_Type = "Gene_Fusion",
    Description = Gene
  ) %>%
  select(colnames(genvar_mut))


genvar_tbl <- bind_rows(genvar_mut, genvar_cop, genvar_fus) %>%
  arrange(TCGA_Cohort, Participant_ID, Sample_ID, 
          Variant_Type, Gene, Description)

stopifnot(
  !any(is.na(genvar_tbl$TCGA_Cohort)),
  !any(is.na(genvar_tbl$Participant_ID)),
  !any(is.na(genvar_tbl$Sample_ID))
)

# OUTPUT -----------------------------------------------------------------------
readr::write_tsv(genvar_tbl, file = gzfile(output_file))
message(paste0("Combined variant data table written to: ", output_file))
message(paste0("Clinical data table written to: ", output_file_clinical))
message(paste0("CNA data table written to: ", output_file_cna))
message(paste0("Fusion data table written to: ", output_file_fusions))
message(paste0("MAF table written to: ", output_file_maf))
# ------------------------------------------------------------------------------



