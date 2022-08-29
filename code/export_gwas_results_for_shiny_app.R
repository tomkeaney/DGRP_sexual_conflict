# Getting some practice data for Becca's app

library(tidyverse)

all_files <- list.files("gwas_data/derived/gwas_results", full.names = T)
all_files <- all_files[grepl("tsv", all_files)]
traits <- str_remove_all(str_remove_all(all_files, "gwas_data/derived/gwas_results/"), "[.]tsv[.]gz")


gwas_hits <- lapply(1:length(all_files), function(i) {
  print(i)
  read_tsv(all_files[i]) %>%
    arrange(P) %>%
    filter(P < 10e-5) %>%
    rename(Variant = SNP, Beta = BETA) %>%
    mutate(Trait = traits[i], .before = 1)
}) %>%
  bind_rows()

write_csv(gwas_hits, "~/Desktop/gwas_hits_for_shiny_app.csv")

manhattan_trait_numbers <- tibble(file_number = 1:length(traits), trait = traits)
write_csv(manhattan_trait_numbers, "~/Desktop/manhattan_plots/trait_numbers.csv")

make_one_manhattan <- function(i, manhattan_trait_numbers, all_files){

  file <- paste("gwas_data/derived/gwas_results/", manhattan_trait_numbers$trait[i], ".tsv.gz", sep = "")

  trait <- str_remove_all(str_remove_all(all_files[i], "gwas_data/derived/gwas_results/"), "[.]tsv[.]gz")
  output_file <- paste("~/Desktop/manhattan_plots/trait_", i, ".png", sep = "")
  print(i)
  print(trait)

  manhattan_data <- read_tsv(all_files[i], col_select = c(SNP, P)) %>%
    mutate(position = str_split(SNP, "_"),
           chr = map_chr(position, ~ .x[1]),
           position = as.numeric(map_chr(position, ~ .x[2]))) %>%
    filter(chr != "4")

  max_pos <- manhattan_data %>%
    group_by(chr) %>%
    summarise(max_pos = max(position), .groups = "drop") %>%
    as.data.frame()
  max_pos$max_pos <- c(0, cumsum(max_pos$max_pos[1:4]))

  manhattan_data <- manhattan_data %>%
    left_join(max_pos, by = "chr") %>%
    mutate(position = position + max_pos)

  plot <- manhattan_data %>%
    ggplot(aes(position, -1 * log10(P), group = chr, fill = chr, stroke = 0.05)) +
    geom_point(size = 0.9, colour="grey20", pch = 21) +
    geom_hline(yintercept = 5, linetype=2) +
    scale_fill_brewer(palette = "Paired", name = "Chromosome") +
    ylab(expression(paste("-", Log[10], " p value"))) +  xlab("Position in the genome") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          legend.position = "top",
          panel.border = element_blank(),
          axis.ticks.x = element_blank())

  ggsave(filename = output_file, plot, height = 9*1.5, width = 14*1.5, units = "cm")
}

lapply(1:nrow(manhattan_trait_numbers), function(i) make_one_manhattan(i, manhattan_trait_numbers, all_files))



# Funky analysis of weighted beta for every SNP

library(tidyverse)

all_files <- list.files("~/Rprojects/DGRP_sexual_conflict/gwas_data/derived/gwas_results", full.names = T)
all_files <- all_files[grepl("tsv", all_files)]
traits <- str_remove_all(str_remove_all(all_files, "gwas_data/derived/gwas_results/"), "[.]tsv[.]gz")
all_snps <- read_tsv(all_files[1]) %>% pull(SNP)
chunked_snps <- split(all_snps, ceiling(seq_along(all_snps)/5000))
n_snps <- length(all_snps)

do_snp_chunk <- function(i){
  print(i)
  focal_snps <- chunked_snps[[i]]
  all_betas <- map_df(1:length(all_files),
                      ~ read_tsv(all_files[.x],
                                 col_types = c("c", "n", "n", "n")) %>%
                        filter(SNP %in% focal_snps) %>%
                        mutate(trait = .x)
  )
  all_betas %>%
    group_by(SNP) %>%
    summarise(Mean_beta = weighted.mean(BETA, w = 1 / SE))
}

weighted_betas <- map_df(1:length(chunked_snps), ~ do_snp_chunk(.x))

