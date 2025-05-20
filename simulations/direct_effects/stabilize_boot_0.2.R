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

source("../../methods.R")

set.seed(321)
nsim = 100 # 100 datasets
bootstrap_reps = 500 # 500 bootstrap repetitions for each

pi <- 0.2
algorithm <- stabilize
overwrite <- TRUE

algorithm_name <- "stabilize"
data_dir <- glue("../datasets_{pi}")
name <- glue("{algorithm_name}_bootstrapped_ATE_{pi}")

if (dir.exists(name)) {
  if (!overwrite) {
    stop(glue("Directory {name} exists and overwrite = {overwrite}"))
  }
} else {
  dir.create(name)
}

results <- matrix(nrow = 0, ncol = 7)
for (i in 1:nsim) {
  sim <- readRDS(file.path(data_dir, glue("sim_{i}.rds")))
  N <- ncol(sim$params$P)
  bootstraped_ATE = bootstrap_wrapper(
    algorithm = algorithm,
    filename = glue("{name}/boot_{i}"),
    reps = bootstrap_reps,
    M = sim$M,
    Tr = sim$Tr,
    N = N,
    reference_P = sim$params$P
  )
  # record bootstrapped ATE and 95% CI
  estimates <- read.csv(glue("{name}/boot_{i}.csv"))
  results_i <- data.frame(
    sim = i,
    mean = apply(estimates, 2, function(col) {mean(col, na.rm = TRUE)}),
    se = apply(estimates, 2, function(col) {sd(col, na.rm = TRUE)}),
    upper95 =  apply(estimates, 2, function(col) {quantile(col, 0.975, na.rm = TRUE)}),
    lower95 = apply(estimates, 2, function(col) {quantile(col, 0.025, na.rm = TRUE)}),
    upper90 =  apply(estimates, 2, function(col) {quantile(col, 0.95, na.rm = TRUE)}),
    lower90 = apply(estimates, 2, function(col) {quantile(col, 0.05, na.rm = TRUE)})
  )
  if (i == 1) {
    results <- results_i
  } else {
    results <- rbind(results, results_i)
  }
  write.csv(results, glue("{name}.csv"), row.names = FALSE)
}
