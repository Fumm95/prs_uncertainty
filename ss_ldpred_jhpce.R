# This script is adapted from https://github.com/privefl/paper-infer/blob/main/code/
# example-with-provided-LD.R

suppressPackageStartupMessages({
  library(bigsnpr)
  library(readr)
  library(tidyverse)
  library(glue)
  library(bigreadr)
  library(matrixStats)
})

################################################################
### Step 1 loading in summary statistics and setting up data
################################################################
# set up parameters

for(chr_num in c(1:22)){
  
  set.seed(246)
  
  ncores <- 4
  n_burn_in <- 1000
  n_chain <- 10
  n_iter <- 500
  report_step <- 5
  pheno <- 2
  
  ncase <- 76192
  ncontrol <- 63082
  ss_file <- paste0('/dcs04/nilanjan/data/mfu/prs_simulation/simulated_betas/bcac_ukbb_P',
                    pheno, 'sumstats.txt')
  
  test_bim_tmp <- '/dcs04/nilanjan/data/mfu/ukbb_bfiles/test_files_small/chrCHROM_test'
  # CHROM will be replaced by chromosome number later
  output_prefix <- paste0('/dcs04/nilanjan/data/mfu/ss_brca/simulated_betas/brca_ukbb_P',
                          pheno, '_out_', chr_num)
  
  ldr = 3/1000
  
  setwd("/fastscratch/myscratch/mfu/")
  
  # load summary statistics
  sumstats <- fread2(ss_file, select = c("chr", "pos", "a0", "a1",
                                         "new_beta", "beta_se", "maf.ss"),
                     col.names = c("chr", "pos", "a0", "a1", "beta", "beta_se", "maf"))
  sumstats <- sumstats %>%
    filter(chr!='X', chr!='Y') %>%
    mutate(chr = as.integer(chr)) %>%
    drop_na(chr) %>% filter(chr==chr_num) # select autosome only
  
  sumstats$n_eff <- 4 / (1 / ncase + 1 / ncontrol)
  
  
  
  ################################################################
  ### Step 2 train PGS 
  ################################################################ 
  
  # create ld matrix 
  tmp <- tempfile(tmpdir = "tmp-data")
  
  obj.bigSNP <- snp_attach(paste0('/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/
                                  chr',chr_num,'_1000.rds'))
  
  G   <- obj.bigSNP$genotypes
  CHR <- obj.bigSNP$map$chromosome
  POS <- obj.bigSNP$map$physical.pos
  
  map <- obj.bigSNP$map[-(2:3)] 
  names(map) <- c("chr", "pos", "a0", "a1")
  
  sumstats$beta <- as.numeric(sumstats$beta)
  info_snp <- snp_match(sumstats, map, strand_flip = T)
  rownames(info_snp) = info_snp$rsid
  
  POS2 <- snp_asGeneticPos(CHR, POS, ncores = 2)
  
  ## indices in info_snp
  ind.chr <- which(info_snp$chr == chr_num)
  df_beta <- info_snp[ind.chr, c("beta", "beta_se", "n_eff")]
  
  nas <- unique(which(is.na(df_beta$beta)), which(is.na(df_beta$beta_se)))
  
  ind.chr <- ind.chr[!(ind.chr %in% nas)]
  
  info_snp_small <- info_snp[-nas,]
  if(length(nas) > 0){df_beta <- df_beta[-nas,]}
  ## indices in G
  ind.chr2 <- info_snp$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(G, ind.col = ind.chr2, ncores = 4,
                   infos.pos = POS2[ind.chr2], size = ldr) # default
  
  corr <- as_SFBM(corr0, tmp, compact = TRUE)
  
  # ldscore to get a starting value 
  (ldsc <- snp_ldsc2(corr0, df_beta))
  h2_est <- abs(ldsc[["h2"]])
  
  # LDPred2 auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.2, length.out = n_chain),
                                 ncores = 5, report_step = report_step, num_iter = n_iter,
                                 burn_in = n_burn_in, allow_jump_sign = FALSE,
                                 shrink_corr = 0.95) 
  
  # rescale post_beta_sample to allelic scale, the default scale of snp_ldpred2_auto output is:
  # if sumstats is on allelic scale, then beta is allelic scale, beta_sample is std scale
  # we need to rescale the effect size to allelic level 
  effect_scales <- with(df_beta, sqrt(n_eff * beta_se^2 + beta^2))  
  for (chain_i in 1:length(multi_auto)){
    multi_auto[[chain_i]]$sample_beta <- sweep(multi_auto[[chain_i]]$sample_beta, 1,
                                               effect_scales, '*')
  }
  
  
  # mini chain quality control
  all_h2 <- sapply(multi_auto, function(auto) auto$h2_est)
  h2 <- median(all_h2)
  keep <- between(all_h2, 0.7 * h2, 1.4 * h2)
  all_p <- sapply(multi_auto, function(auto) auto$p_est)
  p <- median(all_p[keep])
  keep <- keep & between(all_p, 0.5 * p, 2 * p)
  post_beta_sample <- do.call( 'cbind', lapply(multi_auto[keep], function(auto) auto$sample_beta))
  
  colnames(post_beta_sample) <- paste0("SAMPLE_", seq(1, ncol(post_beta_sample)))                                
  
  if(length(nas) > 0){
    post_beta_sample <- bind_cols(
      info_snp_small %>% select(
        chr, pos, a0, a1,
      ),
      as_tibble(as.matrix(post_beta_sample))
    )
  } else{
    post_beta_sample <- bind_cols(
      info_snp %>% select(
        chr, pos, a0, a1,
      ),
      as_tibble(as.matrix(post_beta_sample))
    )
  }
  
  write_tsv(post_beta_sample,  paste0(output_prefix, ".post_beta_sample.tsv.gz"))
  
}