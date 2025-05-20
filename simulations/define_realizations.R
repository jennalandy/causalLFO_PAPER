set.seed(321)
G <- 100
R <- 20

realizations <- list()

for (pi in c(0.2, 0.8)) {
    realizations[[as.character(pi)]] <- list()
    for (r in 1:R) {
        Tr <- sample(c(0,1), G, prob = c(1-pi, pi), rep = TRUE)
        realizations[[as.character(pi)]][[r]] <- Tr
    }
}

saveRDS(realizations, "processed/realizations.rds")