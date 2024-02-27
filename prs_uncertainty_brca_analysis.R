library(data.table)
library(tidyverse)

dat <- fread("/dcs04/nilanjan/data/mfu/ss_brca/test_small/brca_fem_out.post_prs_sample.tsv.gz")

dat <- dat[order(dat$IID),]

prs_means <- dat %>% select(starts_with("SAMPLE")) %>%
  rowMeans()

prs_ses <- dat %>% select(starts_with("SAMPLE")) %>%
  apply(1, sd)

prs_df <- cbind(dat$IID, prs_means, prs_ses)
colnames(prs_df) <- c("ID", "PRS", "SEs")

load("/dcl01/chatterj/data/ameisner/cPRS/outcomePRS20191016.Rdata")

bc_ukbb <- alldata %>% select("subjectID", "sex", "brca_baseline", "brca_incident",
                              "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7",
                              "PC8", "PC9", "PC10")

brca_full <- merge(prs_df, bc_ukbb, by.x = "ID", by.y = "subjectID")

brca_full <- brca_full %>% filter(sex == "Female")

ref_pop <- read.table("/dcs04/nilanjan/data/ydun/PRS_Bridge/ref_ukbb/ldpred_ref/1000.txt",
                      header = F)

all_unrelated <- read.table("/dcl01/chatterj/data/ukbiobank/cleaning_information_data/
                            unrelated_european_ancestry_eid.txt", header = T)

brca_full <- brca_full %>% filter(ID %in% all_unrelated$eid)

brca_full <- brca_full %>% filter(!(ID %in% ref_pop$V1))
brca_inc <- brca_full %>% filter(brca_baseline == 0)

brca_inc$stdprs <- brca_inc$PRS / sd(brca_inc$PRS)

### Mean Relative Risk 

logfit <- glm(brca_incident ~ stdprs + PC1 + PC2 + PC3 +
                PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,
              family = binomial(link = "logit"), data = brca_inc)
summary(logfit)

brca_inc$scprs <- brca_inc$stdprs * (coefficients(logfit)[2])
brca_inc$scvar <- brca_inc$SEs * (coefficients(logfit)[2])

rank_of_mean_prs <- rank(brca_inc$scprs)

### Threshold Exceedence

c <- quantile(brca_inc$scprs, probs = 0.8)
brca_inc$ci_score <- 1 - pnorm(c, mean = brca_inc$scprs, brca_inc$scvar)

rank_of_threshold_exceedence <- rank(brca_inc$ci_score)

### Expected Rank of Relative Risk

prs_mat <- dat %>% select(starts_with("SAMPLE"))

expected_rank_of_prs <- prs_mat %>% apply(2, rank) %>% rowMeans()

### Mean Absolute Risk

brca_inc$abs_risk <- exp(coefficients(logfit)[1] + brca_inc$scprs + 0.5 * brca_inc$scvar^2)

rank_of_abs_risk <- rank(brca_inc$abs_risk)

### Threshold Exceedence of Relative Risk

absfunc <- function(vec){
  return(exp((coefficients(logfit)[1]) + vec * (coefficients(logfit)[2])))
}

absrisk_mat2 <- apply(prs_mat, 2, absfunc)

absrisk_mean2 <- absrisk_mat2 %>% rowMeans()

cemp <- quantile(absrisk_mean2, probs = 0.8)

binmat = absrisk_mat2 > cemp

probvec = rowMeans(binmat)

brca_inc$probemp = probvec

expected_rank_of_absrisk = rank(brca_in$probemp)

### Expected Rank of Absolute Risk 

exp_absrank = apply(absrisk_mat2, 2, rank)

mean_expabsrisk = exp_absrank %>% rowMeans()

expected_rank_of_absolute_risk = rank(absrisk_mean2)

### Comparing Efficacy

brca_inc$exprank_prs = expected_rank_of_prs
brca_inc$exprank_absrisk = expected_rank_of_absolute_risk

brca_full <- brca_inc

ordervec <- rep(0, 6)

ordermat <- brca_full[order(-brca_full$scprs)[c(1:1000)],]
ordervec[1] <- sum(ordermat$brca_incident)

ordermat <- brca_full[order(-brca_full$exprank_prs)[c(1:1000)],]
ordervec[2] <- sum(ordermat$brca_incident)

ordermat <- brca_full[order(-brca_full$ci_score)[c(1:1000)],]
ordervec[3] <- sum(ordermat$brca_incident)

ordermat <- brca_full[order(-brca_full$abs_risk)[c(1:1000)],]
ordervec[4] <- sum(ordermat$brca_incident)

ordermat <- brca_full[order(-brca_full$probemp)[c(1:1000)],]
ordervec[5] <- sum(ordermat$brca_incident)

ordermat <- brca_full[order(-brca_full$exprank_absrisk)[c(1:1000)],]
ordervec[6] <- sum(ordermat$brca_incident)

ordervec