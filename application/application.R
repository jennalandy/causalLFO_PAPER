.libPaths("~/apps/R_4.1.0/")

library(ggplot2, lib = "~/apps/R_4.1.0/")
library(tidyverse, lib = "~/apps/R_4.1.0/")
library(pbapply, lib = "~/apps/R_4.1.0/")
library(BiocGenerics, lib = "~/apps/R_4.1.0/")
library(Biobase, lib = "~/apps/R_4.1.0/")
library(rngtools, lib = "~/apps/R_4.1.0/")
library(registry, lib = "~/apps/R_4.1.0/")
library(NMF, lib = "~/apps/R_4.1.0/")
library(SnowballC, lib = "~/apps/R_4.1.0/")
library(lsa, lib = "~/apps/R_4.1.0/")
library(nnls, lib = "~/apps/R_4.1.0/")
library(glue, lib = "~/apps/R_4.1.0/")

library(RcppHungarian, lib = "~/apps/R_4.1.0/")
library(bayesNMF, lib = "~/apps/R_4.1.0/")

source("../methods.R")
source("../visualization.R")

cosmic <- bayesNMF::get_cosmic()

results_dir <- "application_results"

run <- function(
    M, Tr, rank,
    algorithm, algorithm_name,
    results_dir = results_dir
) {
    bootstrap_res <- bootstrap_wrapper(
        algorithm,
        filename = file.path(results_dir, glue("{algorithm_name}_rank{rank}")),
        M = M,
        Tr = Tr,
        N = rank
    )

    bootstrap_ATEs <- read.csv(file.path(results_dir, glue("{algorithm_name}_rank{rank}.csv")))
    ATE <- bootstrap_ATEs %>% colMeans()
    lower <- apply(bootstrap_ATEs, 2, function(col) quantile(col, 0.025))
    upper <- apply(bootstrap_ATEs, 2, function(col) quantile(col, 0.975))

    aligned_Ps <- readRDS(file.path(results_dir, glue("{algorithm_name}_rank{rank}_aligned_Ps.rds")))
    mean_P <- Reduce(`+`, aligned_Ps) / length(aligned_Ps)
    sim <- pairwise_sim(mean_P, cosmic)
    aligned_sim <- assign_signatures(sim)

    return(data.frame(
      ATE = ATE,
      lower = lower,
      upper = upper,
      sig = colnames(aligned_sim),
      cosine_sim = diag(aligned_sim)
    ))
}

Tr <- read.csv("T_BA_sub45.csv")$BA_T
M <- as.matrix(read.csv("M_BA.csv"))[,paste0("Breast.AdenoCA..", read.csv("T_BA_sub45.csv")$tumor_wgs_icgc_specimen_id)]

methods_list <- list(
  "all_data" = all_data,
  "random_split" = random_split_and_nmf,
  "impute" = impute,
  "stabilize" = stabilize,
  "impute_and_stabilize" = impute_and_stabilize
)
ranks <- 2:5
grid <- expand.grid('method' = names(methods_list), 'rank' = ranks)

# dataset of focus specified in command line argument
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])
method <- grid$method[i]
rank <- grid$rank[i]

print(paste(method, rank))

algorithm <- methods_list[[method]]

this_res <- run(
    M = M,
    Tr = Tr,
    rank = rank,
    algorithm = algorithm,
    algorithm_name = method,
    results_dir = results_dir
)