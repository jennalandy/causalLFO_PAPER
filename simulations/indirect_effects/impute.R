source("../indirect_effects.R")
source("../../methods.R")
realizations <- readRDS("../processed/realizations.rds")

# dataset of focus specified in command line argument
args <- commandArgs(trailingOnly = TRUE)
sim <- as.numeric(args[1])

# using data from pi = 0.2, but doesn't actually use that assigned treatment (uses realizations)
# so could technically use either
data <- readRDS(file.path("../datasets_0.2", glue("sim_{sim}.rds")))
N <- ncol(data$params$P)

algorithm <- impute
dir <- file.path("impute", glue("impute_{sim}"))

# set seed for reproducibility
set.seed(321)
calc_all_IAPOs(
    realizations, 
    algorithm = algorithm, 
    M0 = data$M0, 
    M1 = data$M1, 
    N = N, 
    reference_P = data$params$P, 
    seed = NULL, nrun = 5,
    dir = dir,
    imputed = TRUE
)