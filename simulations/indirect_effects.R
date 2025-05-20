.libPaths("~/apps/R_4.1.0/")
library(glue)
library(stringr)

# This function takes in potential outcomes matrices with columns as subjects
# and a treatment vector, and returns the observed data matrix
combine_mat <- function(M0, M1, Tr) {
  G <- ncol(M1)
  M <- do.call(cbind, lapply(1:G, function(g) {
    if (Tr[g] == 1) {
      return(M1[,g])
    } else {
      return(M0[,g])
    }
  }))
  return(M)
}

# This function creates a directory which may or may not be nested
mk_nested_dir <- function(dir) {
    dir_order <- str_split(dir, '/')[[1]]
    for (i in 1:length(dir_order)) {
        nested_dir = paste(dir_order[1:i], collapse = '/')
        if (!dir.exists(nested_dir)) {
            print(glue("creating {nested_dir}"))
            dir.create(nested_dir)
        }
    }
}

# This function takes in a set of fixed realizations of treatment assignments, 
# an algorithm, and relevant data, and computes estimates of all individual average potential outcomes
# for each subject under each treatment level
calc_all_IAPOs <- function(
    realizations, algorithm, 
    M0, M1, N, 
    reference_P = NULL, 
    seed = NULL, nrun = 5,
    dir = "IAPOs",
    imputed = FALSE,
    ...
) {
    G <- ncol(M0)
    if (!dir.exists(dir)) {
        mk_nested_dir(dir)
    }

    # realizations is a list named by pi values, values are lists of realizations under that pi
    for (pi in names(realizations)) {
        print(glue("pi = {pi}"))

        realizations_pi <- realizations[[as.character(pi)]]
        R <- length(realizations_pi)

        counts_0 <- matrix(0, nrow = N, ncol = G)
        counts_1 <- matrix(0, nrow = N, ncol = G)
        Csum_0 <- matrix(0, nrow = N, ncol = G) # will hold IAPO (sum, not mean) for each individual g (cols) when Tg = 0
        Csum_1 <- matrix(0, nrow = N, ncol = G) # will hold IAPO (sum, not mean) for each individual g (cols) when Tg = 1

        if (imputed) {
            counts_imputed_0 <- matrix(0, nrow = N, ncol = G)
            counts_imputed_1 <- matrix(0, nrow = N, ncol = G)
            Csum_imputed_0 <- matrix(0, nrow = N, ncol = G) # will hold imputed IAPO (sum, not mean) for each individual g (cols) when Tg = 0
            Csum_imputed_1 <- matrix(0, nrow = N, ncol = G) # will hold imputed IAPO (sum, not mean) for each individual g (cols) when Tg = 1
        }
        

        # for each realization, record "observed" outcomes
        if (!isTRUE(all.equal(algorithm, random_split_and_nmf))) {
            # one run under each realization takes care of half of IAPOs
            for (r in 1:R) {
                print(glue("\tr = {r}"))
                Tr <- realizations_pi[[r]]
                M <- combine_mat(M0, M1, Tr)

                est <- algorithm(
                    M, Tr, N, reference_P = reference_P,
                    seed = seed, nrun = nrun,
                    ...
                )
                
                # log the estimated C values as well as the number of times each
                # individual / treatment level combo was recorded

                Csum_0[,Tr == 0] <- Csum_0[,Tr == 0] + est$Chat[,Tr == 0]
                counts_0[,Tr == 0] <- counts_0[,Tr == 0] + 1

                Csum_1[,Tr == 1] <- Csum_1[,Tr == 1] + est$Chat[,Tr == 1]
                counts_1[,Tr == 1] <- counts_1[,Tr == 1] + 1

                if (imputed) {
                    Csum_imputed_0[,Tr == 1] <- Csum_imputed_0[,Tr == 1] + est$Chat_imputed[,Tr == 1] # if Tr = 1, imputed value is for Tr = 0
                    counts_imputed_0[,Tr == 1] <- counts_imputed_0[,Tr == 1] + 1

                    Csum_imputed_1[,Tr == 0] <- Csum_imputed_1[,Tr == 0] + est$Chat_imputed[,Tr == 0] # if Tr = 0, imputed value is for Tr = 1
                    counts_imputed_1[,Tr == 0] <- counts_imputed_1[,Tr == 0] + 1
                }
            }
        } else {
            # need to force_second one at a time if random_split_and_nmf
            # so requires one run per subject/realization combo
            for (g in 1:G) {
                print(glue("\tg = {g}"))
                for (r in 1:R) {
                    print(glue("\t\tr = {r}"))
                    Tr <- realizations_pi[[r]]
                    M <- combine_mat(M0, M1, Tr)

                    est <- algorithm(
                        M, Tr, N, reference_P = reference_P,
                        seed = seed, nrun = nrun,
                        force_second = c(g),
                        ...
                    )
                    
                    if (Tr[g] == 0) {
                        Csum_0[,g] <- Csum_0[,g] + est$Chat[,g]
                        counts_0[,g] <- counts_0[,g] + 1
                        if (imputed) {
                            Csum_imputed_1[,g] <- Csum_imputed_1[,g] + est$Chat_imputed[,g]
                            counts_imputed_1[,g] <- counts_imputed_1[,g] + 1
                        }
                    } else {
                        Csum_1[,g] <- Csum_1[,g] + est$Chat[,g]
                        counts_1[,g] <- counts_1[,g] + 1
                        if (imputed) {
                            Csum_imputed_0[,g] <- Csum_imputed_0[,g] + est$Chat_imputed[,g]
                            counts_imputed_0[,g] <- counts_imputed_0[,g] + 1
                        }
                    }
                }
            }
        }
        # save intermediate results
        write.csv(Csum_0, file.path(dir, glue("Csum_0_pi{pi}.csv")), row.names = FALSE)
        write.csv(Csum_1, file.path(dir, glue("Csum_1_pi{pi}.csv")), row.names = FALSE)
        write.csv(counts_0, file.path(dir, glue("counts_0_pi{pi}.csv")), row.names = FALSE)
        write.csv(counts_1, file.path(dir, glue("counts_1_pi{pi}.csv")), row.names = FALSE)
        if (imputed) {
            write.csv(Csum_imputed_0, file.path(dir, glue("Csum_0_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(Csum_imputed_1, file.path(dir, glue("Csum_1_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(counts_imputed_0, file.path(dir, glue("counts_0_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(counts_imputed_1, file.path(dir, glue("counts_1_imputed_pi{pi}.csv")), row.names = FALSE)
        }
        
        # for each each realization, record "unobserved" outcomes
        # these require flipping each treatment one at a time
        for (g in 1:G) {
            print(glue("\tg = {g}"))
            for (r in 1:R) {
                print(glue("\t\tr = {r}"))
                Tr <- realizations_pi[[r]]

                # flip treatment of sample g
                Tr[g] <- 1 - Tr[g]
                M <- combine_mat(M0, M1, Tr)

                est <- algorithm(
                    M, Tr, N, reference_P = reference_P,
                    seed = seed, nrun = nrun,
                    force_second = c(g),
                    ...
                )

                if (Tr[g] == 0) {
                    Csum_0[,g] <- Csum_0[,g] + est$Chat[,g]
                    counts_0[,g] <- counts_0[,g] + 1
                    if (imputed) { # if Tr = 0, imputed value is for Tr = 1
                        Csum_imputed_1[,g] <- Csum_imputed_1[,g] + est$Chat_imputed[,g]
                        counts_imputed_1[,g] <- counts_imputed_1[,g] + 1
                    }
                } else {
                    Csum_1[,g] <- Csum_1[,g] + est$Chat[,g]
                    counts_1[,g] <- counts_1[,g] + 1
                    if (imputed) { # if Tr = 1, imputed value is for Tr = 0
                        Csum_imputed_0[,g] <- Csum_imputed_0[,g] + est$Chat_imputed[,g]
                        counts_imputed_0[,g] <- counts_imputed_0[,g] + 1
                    }
                }
            }
        }

        # save final results
        write.csv(Csum_0, file.path(dir, glue("Csum_0_pi{pi}.csv")), row.names = FALSE)
        write.csv(Csum_1, file.path(dir, glue("Csum_1_pi{pi}.csv")), row.names = FALSE)
        write.csv(counts_0, file.path(dir, glue("counts_0_pi{pi}.csv")), row.names = FALSE)
        write.csv(counts_1, file.path(dir, glue("counts_1_pi{pi}.csv")), row.names = FALSE)
        if (imputed) {
            write.csv(Csum_imputed_0, file.path(dir, glue("Csum_0_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(Csum_imputed_1, file.path(dir, glue("Csum_1_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(counts_imputed_0, file.path(dir, glue("counts_0_imputed_pi{pi}.csv")), row.names = FALSE)
            write.csv(counts_imputed_1, file.path(dir, glue("counts_1_imputed_pi{pi}.csv")), row.names = FALSE)
        }
    }
}