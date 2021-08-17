## ----include=FALSE, eval=FALSE-------------------------------------------------------------------------------------------------------------------
## Note from Luke to Luke: To run mashr on Spartan, use `knitr::purl(input = "analysis/run_mashr.Rmd", output = "code/run_mashr.R")` to generate an R script from this R Markdown document, then type this on Spartan
## 
## cd /data/projects/punim0243/DGRP_sexual_conflict
## module load r/4.0.0
## Rscript code/run_mashr.R


## ------------------------------------------------------------------------------------------------------------------------------------------------
#setwd("/data/projects/punim0243/DGRP_mashr")

library(tidyverse)
library(ashr) 
library(mashr) 
library(glue)
library(rslurm)


files <- list.files("gwas_data/derived/gwas_results", pattern = "tsv.gz", full.names = TRUE)
# fitness_files <- files[grepl("fitness.early", files)] # for testing!!
# files <- c(files[1:5], fitness_files)                 # for testing!!
traits <- str_split(files, "/") %>% map_chr(~ .x[4]) %>% str_remove_all(".tsv.gz")
nonfitness_traits <- traits[!(traits %in% c("fitness.early.f", "fitness.early.m"))]
  
# 
# all_dat <- lapply(1:length(files), function(i){
#   trait <- traits[i]
#   dat <- read_tsv(files[i]) %>% select(SNP, BETA, SE)
#   names(dat)[2] <- paste("beta", trait, sep = "_")
#   names(dat)[3] <- paste("se", trait, sep = "_")
#   dat
# })
# names(all_dat) <- traits
# 
# all_dat <- lapply(all_dat, function(x) head(x, 10000)) # for testing! 


## ------------------------------------------------------------------------------------------------------------------------------------------------
run_mashr <- function(nonfitness_trait, overwrite = FALSE){
  
  setwd("/data/projects/punim0243/DGRP_sexual_conflict")
  print(nonfitness_trait)
  
  # Set up the output file name. If it already exists and overwrite = FALSE, just quit.
  # nonfitness_trait_neat <- str_replace_all(nonfitness_trait, "[.]", "_")
  output_file <- glue("gwas_data/derived/mashr_results/{nonfitness_trait}.rds")
  if(!overwrite & file.exists(output_file)) return(NULL)
  
  # Input file names (GWAS stats for 2 fitness traits and a another phenotype)
  files <- c(
    "gwas_data/derived/gwas_results/fitness.early.f.tsv.gz",
    "gwas_data/derived/gwas_results/fitness.early.m.tsv.gz",
    paste("gwas_data/derived/gwas_results/", nonfitness_trait, ".tsv.gz", sep = ""))
  traits <- c("fitness.early.f", "fitness.early.m", "nonfitness_trait")
  
  # Load the GWAS data for these 3 things
  all_dat <- lapply(1:length(files), function(i){
    trait <- traits[i]
    dat <- read_tsv(files[i]) %>% select(SNP, BETA, SE)
    names(dat)[2] <- paste("beta", trait, sep = "_")
    names(dat)[3] <- paste("se", trait, sep = "_")
    dat
  })
  names(all_dat) <- traits

  # all_dat <- lapply(all_dat, function(x) head(x, 10000)) # for testing! 
  
  ##### 1. Data set up
  # Define functions to get the focal data, and set up as betas and SE for mashr
  make_mashr_data_one_phenotype <- function(nonfitness_trait){
    left_join(all_dat[[1]], all_dat[[2]], by = "SNP") %>% # 1=female fitness, 2=male fitness, 3=phenotype
      left_join(all_dat[[3]], by = "SNP") %>% 
      select(SNP, starts_with("beta"), starts_with("SE"))
  }
  
  mashr_setup <- function(beta_and_se){
    betas <- beta_and_se %>% select(starts_with("beta")) %>% as.matrix()
    SEs <- beta_and_se %>% select(starts_with("SE")) %>% as.matrix()
    rownames(betas) <- beta_and_se$SNP
    rownames(SEs) <- beta_and_se$SNP
    mash_set_data(betas, SEs)
  }
  
  beta_and_se <- make_mashr_data_one_phenotype(nonfitness_trait)
  mash_data <- mashr_setup(beta_and_se)
  
  ##### 2. Obtain data-driven covariance matrices for mashr (see mashr vignette)
  # This first finds strong signals in the complete SNP effect size data, and use
  # to estimate the underlying covariance matrices for the SNPs' effect sizes
  m.1by1 <- mash_1by1(mash_data) 
  strong <- get_significant_results(m.1by1, thresh = 0.2)   
  U.pca <- cov_pca(mash_data, npc = 3, subset = strong)
  U <- cov_ed(mash_data, U.pca, subset = strong)
  
  ##### 3. Run mashr and save the results
  # Now that the data and cov matrices are set up, we can run mashr (takes an hour or so?)
  mash_ED <- mash(data = mash_data, Ulist = U)
  list(SNPs = beta_and_se %>% pull(SNP), mash_ED) %>% saveRDS(output_file)
}

# mashr_results <- run_mashr("aggression.m")

slurm_parameters_df <- data.frame(nonfitness_trait = nonfitness_traits)
# slurm_parameters_df <- head(slurm_parameters_df) # for testing! 
sopt1 <- list(time = '6:00:00')

sjob <- slurm_apply(f = run_mashr, params = slurm_parameters_df, 
                    jobname = 'DGRP_mashr', slurm_options = sopt1,
                    nodes = nrow(slurm_parameters_df), cpus_per_node = 1, submit = TRUE)

# sort(get_estimated_pi(mash_ED)) %>%
#   enframe() %>% arrange(value) %>% 
#   filter(value > 0.01) %>% 
#   mutate(name = factor(name, unique(name))) %>% 
#   ggplot(aes(name, 100 * value)) + 
#   geom_bar(stat = "identity", colour = "grey10",  size = 0.3) + 
#   coord_flip() 

