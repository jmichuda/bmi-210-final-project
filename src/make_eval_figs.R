library(ggpubr)
library(tidyverse)

# PARAMS ----------------------------------------------------------------------
inference_dir  <- "inference"
txrec_freq_plot <- "figures/txrec_freq_by_cantype.pdf"

# LOAD DATA -------------------------------------------------------------------
recfiles <- list.files(inference_dir, full.names = TRUE, 
                       pattern = "*tcga_samples_[a-z]+\\.csv")

recdata <- purrr::map(recfiles, readr::read_csv, 
                      show_col_types = FALSE,
                      col_types = list(level = col_character())) %>%
  bind_rows() %>%
  unique() %>%
  mutate(
    therapy = stringr::str_replace(therapy, "^oncokb_civic\\.", ""),
    type = ifelse(level %in% c("R1", "R2"), "RES", "SEN"),
    #rec_id = paste(disease, type, gene, variant, therapy, sep = "_"),
    rec_id = paste(paste0("%", disease, "%"), # to be trimmed below
                   gene, "_", variant, 
                   paste(" (", ifelse(type == "SEN", "+", "-"), ")", sep = ""), 
                   therapy)
  ) %>%
  filter(level != "4") %>% # exclude pre-clinical associations
  select(-type) %>%
  select(rec_id, everything()) %>%
  arrange(rec_id, patient, source)

# TABULATE THERAPEUTIC RECOMMENDATION FREQUENCIES -----------------------------
rec_coverage <- full_join(
  recdata %>% 
    filter(source == "oncokb") %>% 
    count(disease, rec_id, name = "n_oncokb"),
  recdata %>% 
    filter(source == "civic") %>% 
    count(disease, rec_id, name = "n_civic"),
  by = c("disease", "rec_id")
) %>%
  replace_na(list(n_oncokb = 0, n_civic = 0)) %>%
  mutate(diff_okb_civ = n_oncokb - n_civic)

stopifnot(!any(duplicated(rec_coverage$rec_id)))

# PLOT UTILS  -----------------------------------------------------------------
prepare_plot_data <- function(rec_cov_tbl, cantype, n = 10) {
  source_name = c(oncokb = "OncoKB", civic = "CIViC")
  
  tb <- rec_cov_tbl %>% 
    filter(disease == cantype) %>% 
    arrange(desc(n_oncokb), desc(n_civic))
  tb <- tb[1:n, ] %>%
    mutate(rec_id = stringr::str_replace(rec_id, "^%[A-Z]+% ", "")) #rm disease
  ordered_recs <- rev(tb$rec_id)
  
  tb %>%
    select(-diff_okb_civ) %>%
    rename(oncokb = n_oncokb, civic = n_civic) %>%
    pivot_longer(cols = c("oncokb", "civic"), 
                 names_to = "source", values_to = "n") %>%
    mutate(
      rec_id = factor(rec_id, levels = ordered_recs),
      source = source_name[source],
      source = factor(source, levels = c("OncoKB", "CIViC")),
    ) %>%
    select(disease, everything())
}

get_plot_obj <- function(rec_cov_tbl, cantype, n = 10) {
  tb <- prepare_plot_data(rec_coverage, cantype, n)
  ggplot(tb, aes(x = rec_id, y = n, fill = source)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    xlab("") + ylab("# Patients") + 
    ggtitle(unique(tb$disease)) +
    theme(plot.title = element_text(hjust = 0.5))
}

# MAKE FIGURES ----------------------------------------------------------------
p_luad <- get_plot_obj(rec_coverage, "LUAD", n = 20) 

p_lusc <- get_plot_obj(rec_coverage, "LUSC", n = 20) 

p_coad <- get_plot_obj(rec_coverage, "COAD", n = 20)

p_brca <- get_plot_obj(rec_coverage, "BRCA", n = 20)

p_blca <- get_plot_obj(rec_coverage, "BLCA", n = 20)

p_paad <- get_plot_obj(rec_coverage, "PAAD", n = 20)

p <- ggpubr::ggarrange(p_brca, p_luad, 
                       p_coad, p_lusc,
                       p_blca, p_paad,
                       nrow=3, ncol=2, align = "v", 
                       common.legend = TRUE, legend="bottom")

ggsave(plot = p, filename = txrec_freq_plot, width = 11, height = 11)
