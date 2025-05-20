.libPaths("~/apps/R_4.1.0/")
library(NMF)
library(bayesNMF)
library(glue)

oracle <- function(M, Tr, N, reference_P, sim, nnlm_likelihood = NULL) {
  psi_hat_C = rowMeans_wrapper(sim$C1 - sim$C0)
  return(list(
    psi_hat_C = psi_hat_C,
    Phat = reference_P
  ))
}

observed_C <- function(M, Tr, N, reference_P, sim, nnlm_likelihood = NULL) {
  psi_hat_C = rowMeans_wrapper(sim$C1[,Tr == 1]) - rowMeans_wrapper(sim$C0[,Tr == 0])
  return(list(
    psi_hat_C = psi_hat_C,
    Phat = reference_P
  ))
}

all_data <- function(M, Tr, N, reference_P = NULL, seed = NULL, nrun = 5, force_second = NULL, nnlm_likelihood = NULL, bayesian = FALSE, filename = NULL) {
  # NMF on all data
  # this function rescales and aligns to reference
  nmf_res <- nmf_wrapper(
    M, rank = N, nrun = nrun, seed = seed,
    method = 'brunet', reference_P = reference_P,
    bayesian = bayesian,
    filename = paste0(filename, "_nmf")
  )

  # Difference of means estimate of ATE on all data
  psi_hat_C <- rowMeans_wrapper(nmf_res$C[,Tr==1]) -
               rowMeans_wrapper(nmf_res$C[,Tr==0])

  return(list(
    psi_hat_C = psi_hat_C,
    sim_mat = nmf_res$reassigned,
    Chat = nmf_res$C,
    Phat = nmf_res$P
  ))
}

random_split_and_nmf <- function(M, Tr, N, reference_P = NULL, prop = 0.5, seed = NULL, nrun = 5, force_second = c(), nnlm_likelihood = 'poisson', bayesian = FALSE, filename = NULL) {
  split_idx = split_dat(M, Tr, prop, seed = seed, force_second = force_second)

  # need at least one treated and one untreated in the second half
  while (mean(Tr[-split_idx]) == 1 | mean(Tr[-split_idx]) == 0) {
    print("trying split again")
    split_idx = split_dat(M, Tr, prop, seed = seed, force_second = force_second)
  }

  # NMF on first part of data
  # this function rescales and aligns to reference
  nmf_res <- nmf_wrapper(
    M[,split_idx], rank = N, nrun = nrun, seed = seed,
    method = 'brunet', reference_P = reference_P, 
    bayesian = bayesian,
    filename = paste0(filename, "_nmf")
  )

  # Difference of means estimate of ATE on second part of data
  # first NNLS with fixed P to estimate C
  Chat <- nnlm_wrapper(M[,-split_idx], fixed_P = nmf_res$P, likelihood = nnlm_likelihood, bayesian = bayesian, filename = paste0(filename, "_nnlm"))
  psi_hat_C <- rowMeans_wrapper(Chat[,Tr[-split_idx]==1]) -
               rowMeans_wrapper(Chat[,Tr[-split_idx]==0])

  Chat_expanded <- matrix(NA, nrow = N, ncol = ncol(M))
  Chat_expanded[,-split_idx] <- Chat

  return(list(
    psi_hat_C = psi_hat_C,
    sim_mat = nmf_res$reassigned,
    Chat = Chat_expanded,
    split_idx = split_idx,
    Phat = nmf_res$P
  ))
}

calc_psi_hat_M <- function(M, Tr, psi_hat_M_method) {
  if (psi_hat_M_method == "diffMeans") {
    # Difference of means estimate of ATE of M
    psi_hat_M <- rowMeans_wrapper(M[,Tr==1]) -
      rowMeans_wrapper(M[,Tr==0])
  } else if (psi_hat_M_method == "diffMedians") {
    # Difference of medians estimate of ATE of M
    psi_hat_M <- rowMedians_wrapper(M[,Tr==1]) -
      rowMedians_wrapper(M[,Tr==0])
  } else {
    stop(paste("psi_hat_M_method {psi_hat_M_method} not implemented"))
  }
}

impute_and_stabilize <- function(
  M, Tr, N, reference_P = NULL,
  seed = NULL, nrun = 5,
  variance_stabilization = "sqrt",
  psi_hat_M_method = "diffMeans",
  nmf_on = "untreated",
  psi_hat_C_method = "meanITE",
  true_psi_M = NULL,
  force_second = NULL,
  bias_correction = "theoretical",
  nnlm_likelihood = "poisson",
  bayesian = FALSE,
  filename = NULL
) {
  G <- ncol(M)

  if (variance_stabilization == "identity") {

    if (!is.null(true_psi_M)) {
      psi_hat_M = true_psi_M
    } else {
      psi_hat_M = calc_psi_hat_M(M, Tr, psi_hat_M_method)
    }

    # Impute potential outcomes
    M1 <- do.call(cbind, lapply(1:G, function(g) {
      M[,g] + psi_hat_M*(1-Tr[g])
    }))
    M0 <- do.call(cbind, lapply(1:G, function(g) {
      M[,g] - psi_hat_M*Tr[g]
    }))

  } else if (variance_stabilization == "sqrt" | variance_stabilization == "sqrt_bias_corrected") {
    # Perform variance stabilizing
    M_star <- sqrt(M)

    if (!is.null(true_psi_M)) {
      psi_hat_M_star = true_psi_M
      warning(glue("with variance_stabilization sqrt, make sure provided true_psi_M is on correct scale (difference of sqrt potential outcomes)"))
    } else {
      psi_hat_M_star = calc_psi_hat_M(M_star, Tr, psi_hat_M_method)
    }


    # Impute potential outcomes
    M1_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] + psi_hat_M_star*(1-Tr[g])
    }))
    M0_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] - psi_hat_M_star*Tr[g]
    }))

    # Back-transform potential outcomes
    M1 <- M1_star**2
    M0 <- M0_star**2

    # Bias-correction imputed values only, M1[,Tr == 0] and M0[,Tr == 1]
    if (variance_stabilization == "sqrt_bias_corrected") {
      if (bias_correction == 'var') {
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + apply(M1_star, 1, var) # R adds vector to each column
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + apply(M0_star, 1, var)
      } else if (bias_correction == 'observable_var') {
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + apply(M1_star[,Tr == 1], 1, var)
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + apply(M0_star[,Tr == 0], 1, var)
      } else if (bias_correction == 'theoretical') {
        theoretical_correction = 0.25*(1 + 1/sum(Tr == 0) + 1/sum(Tr == 1))
        print(glue("theoretical bias correction {theoretical_correction}"))
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + theoretical_correction
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + theoretical_correction
      }
    }    
  } else if (variance_stabilization == "anscombe") {
    # Perform variance stabilizing
    M_star <- 2*sqrt(M + 3/8)

    if (!is.null(true_psi_M)) {
      psi_hat_M_star = true_psi_M
      warning(glue("with variance_stabilization sqrt, make sure provided true_psi_M is on correct scale (difference of sqrt potential outcomes)"))
    } else {
      psi_hat_M_star = calc_psi_hat_M(M_star, Tr, psi_hat_M_method)
    }


    # Impute potential outcomes
    M1_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] + psi_hat_M_star*(1-Tr[g])
    }))
    M0_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] - psi_hat_M_star*Tr[g]
    }))

    # Back-transform potential outcomes
    M1 <- (2*M1_star)**2 - 3/8
    M0 <- (2*M0_star)**2 - 3/8
    if (sum(M1 < 0) > 0) {
      print(glue("{sum(M1 < 0)} M1 < 0"))
      M1[M1 < 0] <- 0
    }
    if (sum(M0 < 0) > 0) {
      print(glue("{sum(M0 < 0)} M0 < 0"))
      M0[M0 < 0] <- 0
    }
  } else if (variance_stabilization == "poisson_regression") {
    # fit Poisson regression Mg ~ Poisson(lambda(Tg))
    K <- nrow(M)
    lambda_0 <- matrix(NA, nrow = K, ncol = G)
    lambda_1 <- matrix(NA, nrow = K, ncol = G)

    burden <- colSums(M)
    
    M_centered <- scale(t(M), center = TRUE, scale = FALSE)
    pca <- prcomp(M_centered, rank. = 5)
    PCs <- pca$x

    for (k in 1:K) {
      fit <- glm(M[k,] ~ Tr + burden + PCs, family = poisson(link = "log"))
      
      # Predict lambda under T = 0
      Tr0 <- rep(0, G)
      lambda_0[k, ] <- predict(fit, newdata = data.frame(Tr = Tr0, burden = burden, PCs = PCs), type = "response")
      
      # Predict lambda under T = 1
      Tr1 <- rep(1, G)
      lambda_1[k, ] <- predict(fit, newdata = data.frame(Tr = Tr1, burden = burden, PCs = PCs), type = "response")
    }
   
    M0 <- do.call(cbind, lapply(1:G, function(g) {
      rpois(K, lambda = lambda_0[,g])
    }))
    M1 <- do.call(cbind, lapply(1:G, function(g) {
      rpois(K, lambda = lambda_1[,g])
    }))
  } else {
    stop(glue("variance_stabilization {variance_stabilization} not implemented"))
  }

  if (nmf_on == "all") {
    # NMF on all potential outcomes
    M_concat <- cbind(M0, M1)
    nmf_res <- nmf_wrapper(
      M_concat, rank = N, nrun = nrun, seed = seed,
      method = 'brunet', reference_P = reference_P,
      bayesian = bayesian,
      filename = paste0(filename, "_nmf")
    )
    Chat1 <- nmf_res$C[,(G+1):(2*G)]
    Chat0 <- nmf_res$C[,1:G]
  } else if (nmf_on == "untreated") {
    nmf_res <- nmf_wrapper(
      M0, rank = N, nrun = nrun, seed = seed,
      method = 'brunet', reference_P = reference_P,
      bayesian = bayesian,
      filename = paste0(filename, "_nmf")
    )
    Chat0 <- nmf_res$C
    Chat1 <- nnlm_wrapper(M1, fixed_P = nmf_res$P, likelihood = nnlm_likelihood, bayesian = bayesian, filename = paste0(filename, "_nnlm"))
  } else if (nmf_on == "treated") {
    nmf_res <- nmf_wrapper(
      M1, rank = N, nrun = nrun, seed = seed,
      method = 'brunet', reference_P = reference_P,
      bayesian = bayesian,
      filename = paste0(filename, "_nmf")
    )
    Chat1 <- nmf_res$C
    Chat0 <- nnlm_wrapper(M0, fixed_P = nmf_res$P, likelihood = nnlm_likelihood, bayesian = bayesian, filename = paste0(filename, "_nnlm"))
  } else {
    stop(glue("nmf_on {nmf_on} not implemented"))
  }

  # Difference of means estimate of ATE on all data
  if (psi_hat_C_method == "meanITE") {
    # use all estimated potential outcomes of C
    # same sample size so same as difference of means
    # perhaps less noisy because cancels out w/in sample outliers
    psi_hat_C <- rowMeans_wrapper(Chat1 - Chat0)
  } else if (psi_hat_C_method == "diffObsMeans") {
    # only use estimated Chat for the observed treatment level
    # perhaps less noisy as less dependent on imputation
    psi_hat_C <- rowMeans_wrapper(Chat1[,Tr==1]) -
                 rowMeans_wrapper(Chat0[,Tr==0])
  } else if (psi_hat_C_method == "decompose") {
    print("decompose")
    psi_hat_M <- rowMeans_wrapper(M1 - M0)
    psi_hat_C <- decompose(psi_M = psi_hat_M, P = nmf_res$P)
  } else {
    stop(glue("psi_hat_C_method {psi_hat_C_method} not implemented"))
  }

  return(list(
    psi_hat_C = psi_hat_C,
    sim_mat = nmf_res$reassigned,
    Chat = combine_mat(Chat0, Chat1, Tr),
    Chat_imputed = combine_mat(Chat0, Chat1, 1-Tr),
    Phat = nmf_res$P
  ))
}


impute <- function(
  M, Tr, N, reference_P = NULL,
  seed = NULL, nrun = 5,
  variance_stabilization = "sqrt",
  psi_hat_M_method = "diffMeans",
  psi_hat_C_method = "meanITE",
  true_psi_M = NULL,
  force_second = NULL,
  bias_correction = "theoretical",
  nnlm_likelihood = "poisson",
  bayesian = FALSE,
  filename = NULL
) {
  G <- ncol(M)

  if (variance_stabilization == "identity") {

    if (!is.null(true_psi_M)) {
      psi_hat_M = true_psi_M
    } else {
      psi_hat_M = calc_psi_hat_M(M, Tr, psi_hat_M_method)
    }

    # Impute potential outcomes
    M1 <- do.call(cbind, lapply(1:G, function(g) {
      M[,g] + psi_hat_M*(1-Tr[g])
    }))
    M0 <- do.call(cbind, lapply(1:G, function(g) {
      M[,g] - psi_hat_M*Tr[g]
    }))

  } else if (variance_stabilization == "sqrt" | variance_stabilization == "sqrt_bias_corrected") {
    # Perform variance stabilizing
    M_star <- sqrt(M)

    if (!is.null(true_psi_M)) {
      psi_hat_M_star = true_psi_M
      warning(glue("with variance_stabilization sqrt, make sure provided true_psi_M is on correct scale (difference of sqrt potential outcomes)"))
    } else {
      psi_hat_M_star = calc_psi_hat_M(M_star, Tr, psi_hat_M_method)
    }


    # Impute potential outcomes
    M1_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] + psi_hat_M_star*(1-Tr[g])
    }))
    M0_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] - psi_hat_M_star*Tr[g]
    }))

    # Back-transform potential outcomes
    M1 <- M1_star**2
    M0 <- M0_star**2

    # Bias-correction imputed values only, M1[,Tr == 0] and M0[,Tr == 1]
    if (variance_stabilization == "sqrt_bias_corrected") {
      if (bias_correction == 'var') {
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + apply(M1_star, 1, var) # R adds vector to each column
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + apply(M0_star, 1, var)
      } else if (bias_correction == 'observable_var') {
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + apply(M1_star[,Tr == 1], 1, var)
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + apply(M0_star[,Tr == 0], 1, var)
      } else if (bias_correction == 'theoretical') {
        theoretical_correction = 0.25*(1 + 1/sum(Tr == 0) + 1/sum(Tr == 1))
        print(glue("theoretical bias correction {theoretical_correction}"))
        M1[,Tr == 0] <- M1_star[,Tr == 0]**2 + theoretical_correction
        M0[,Tr == 1] <- M0_star[,Tr == 1]**2 + theoretical_correction
      }
    }    
  } else if (variance_stabilization == "anscombe") {
    # Perform variance stabilizing
    M_star <- 2*sqrt(M + 3/8)

    if (!is.null(true_psi_M)) {
      psi_hat_M_star = true_psi_M
      warning(glue("with variance_stabilization sqrt, make sure provided true_psi_M is on correct scale (difference of sqrt potential outcomes)"))
    } else {
      psi_hat_M_star = calc_psi_hat_M(M_star, Tr, psi_hat_M_method)
    }


    # Impute potential outcomes
    M1_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] + psi_hat_M_star*(1-Tr[g])
    }))
    M0_star <- do.call(cbind, lapply(1:G, function(g) {
      M_star[,g] - psi_hat_M_star*Tr[g]
    }))

    # Back-transform potential outcomes
    M1 <- (2*M1_star)**2 - 3/8
    M0 <- (2*M0_star)**2 - 3/8
    if (sum(M1 < 0) > 0) {
      print(glue("{sum(M1 < 0)} M1 < 0"))
      M1[M1 < 0] <- 0
    }
    if (sum(M0 < 0) > 0) {
      print(glue("{sum(M0 < 0)} M0 < 0"))
      M0[M0 < 0] <- 0
    }
  } else if (variance_stabilization == "poisson_regression") {
    # fit Poisson regression Mg ~ Poisson(lambda(Tg))
    K <- nrow(M)
    lambda_0 <- matrix(NA, nrow = K, ncol = G)
    lambda_1 <- matrix(NA, nrow = K, ncol = G)

    burden <- colSums(M)
    
    M_centered <- scale(t(M), center = TRUE, scale = FALSE)
    pca <- prcomp(M_centered, rank. = 5)
    PCs <- pca$x

    for (k in 1:K) {
      fit <- glm(M[k,] ~ Tr + burden + PCs, family = poisson(link = "log"))
      
      # Predict lambda under T = 0
      Tr0 <- rep(0, G)
      lambda_0[k, ] <- predict(fit, newdata = data.frame(Tr = Tr0, burden = burden, PCs = PCs), type = "response")
      
      # Predict lambda under T = 1
      Tr1 <- rep(1, G)
      lambda_1[k, ] <- predict(fit, newdata = data.frame(Tr = Tr1, burden = burden, PCs = PCs), type = "response")
    }
   
    M0 <- do.call(cbind, lapply(1:G, function(g) {
      rpois(K, lambda = lambda_0[,g])
    }))
    M1 <- do.call(cbind, lapply(1:G, function(g) {
      rpois(K, lambda = lambda_1[,g])
    }))
  } else {
    stop(glue("variance_stabilization {variance_stabilization} not implemented"))
  }

  nmf_res <- nmf_wrapper(
    M, rank = N, nrun = nrun, seed = seed,
    method = 'brunet', reference_P = reference_P,
    bayesian = bayesian, filename = paste0(filename, "_nmf")
  )
  Chat0 <- matrix(nrow = nrow(nmf_res$C), ncol = ncol(nmf_res$C))
  Chat1 <- matrix(nrow = nrow(nmf_res$C), ncol = ncol(nmf_res$C))

  Chat0[,Tr == 0] <- nmf_res$C[,Tr == 0]
  Chat1[,Tr == 1] <- nmf_res$C[,Tr == 1]

  MnotT = combine_mat(M0, M1, 1-Tr)
  ChatnotT = nnlm_wrapper(MnotT, fixed_P = nmf_res$P, likelihood = nnlm_likelihood, bayesian = bayesian, filename = paste0(filename, "_nnlm"))
  Chat0[,Tr == 1] <- ChatnotT[,Tr == 1]
  Chat1[,Tr == 0] <- ChatnotT[,Tr == 0]

  # Difference of means estimate of ATE on all data
  if (psi_hat_C_method == "meanITE") {
    # use all estimated potential outcomes of C
    # same sample size so same as difference of means
    # perhaps less noisy because cancels out w/in sample outliers
    psi_hat_C <- rowMeans_wrapper(Chat1 - Chat0)
  } else if (psi_hat_C_method == "diffObsMeans") {
    # only use estimated Chat for the observed treatment level
    # perhaps less noisy as less dependent on imputation
    psi_hat_C <- rowMeans_wrapper(Chat1[,Tr==1]) -
                 rowMeans_wrapper(Chat0[,Tr==0])
  } else {
    stop(glue("psi_hat_C_method {psi_hat_C_method} not implemented"))
  }

  return(list(
    psi_hat_C = psi_hat_C,
    sim_mat = nmf_res$reassigned,
    Chat = combine_mat(Chat0, Chat1, Tr),
    Chat_imputed = combine_mat(Chat0, Chat1, 1-Tr),
    Phat = nmf_res$P
  ))
}

stabilize <- function(
  M, Tr, N, reference_P = NULL,
  seed = NULL, nrun = 5,
  force_second = NULL,
  nnlm_likelihood = "poisson",
  bayesian = FALSE,
  filename = NULL
) {
  G <- ncol(M)

  nmf_0 <- nmf_wrapper(
    M[, Tr == 0], rank = N, nrun = nrun, seed = seed,
    method = 'brunet', reference_P = reference_P,
    bayesian = bayesian, filename = paste0(filename, "_nmf")
  )
  M1 = M[,Tr == 1]
  if (sum(Tr) == 1) {
    M1 = matrix(M1, ncol = 1)
  }
  nnlm_1 <- nnlm_wrapper(
    M1, fixed_P = nmf_0$P, likelihood = nnlm_likelihood,
    bayesian = bayesian, filename = paste0(filename, "_nnlm")
  )

  Chat <- matrix(nrow = N, ncol = G)
  Chat[,Tr == 0] <- nmf_0$C
  Chat[,Tr == 1] <- nnlm_1

  psi_hat_C <- rowMeans_wrapper(Chat[,Tr==1]) -
               rowMeans_wrapper(Chat[,Tr==0])

  return(list(
    psi_hat_C = psi_hat_C,
    sim_mat = nmf_0$reassigned,
    Chat = Chat,
    Phat = nmf_0$P
  ))
}


bootstrap_data <- function(M, Tr) {
  G = ncol(M)
  idx = sample(1:G, G, replace = TRUE)

  return(list(
    M = M[,idx],
    Tr = Tr[idx]
  ))
}

bootstrap_sim <- function(sim) {
  G = ncol(sim$M)
  idx = sample(1:G, G, replace = TRUE)

  return(list(
    M = sim$M[,idx],
    C = sim$C[,idx],
    Tr = sim$Tr[idx],
    M0 = sim$M0[,idx],
    M1 = sim$M1[,idx],
    C0 = sim$C0[,idx],
    C1 = sim$C1[,idx]
  ))
}

bootstrap_wrapper <- function(
  algorithm, filename, M, Tr, N,
  reference_P = NULL,
  reps = 500,
  nnlm_likelihood = 'poisson',
  sim = NULL,
  ...
) {
  if (!is.null(sim)) {
    print("sim provided")
  }
  external_reference <- !is.null(reference_P)
  estimates <- matrix(nrow = reps, ncol = N)
  aligned_Ps <- list()
  for (rep in 1:reps) {
    if (!is.null(sim)) {
      boot_data <- bootstrap_sim(sim)
    } else {
      boot_data <- bootstrap_data(M, Tr)
    }
    while(mean(boot_data$Tr) == 0 | mean(boot_data$Tr) == 1) {
      print("trying bootstrap again")
      if (!is.null(sim)) {
        boot_data <- bootstrap_sim(sim)
      } else {
        boot_data <- bootstrap_data(M, Tr)
      }
    }

    print("boot_data$Tr")
    print(boot_data$Tr)

    # Call the algorithm extra arguments `...`
    if (!is.null(sim)) {
      est <- algorithm(
        M = boot_data$M,
        Tr = boot_data$Tr,
        N = N,
        reference_P = reference_P,
        nnlm_likelihood = nnlm_likelihood,
        sim = boot_data,
        ... # e.g., psi_hat_C_method, FPT, etc.
      )
    } else {
      est <- algorithm(
        M = boot_data$M,
        Tr = boot_data$Tr,
        N = N,
        reference_P = reference_P,
        nnlm_likelihood = nnlm_likelihood,
        ... # e.g., psi_hat_C_method, FPT, etc.
      )
    }
    estimates[rep,] <- est$psi_hat_C
    aligned_Ps[[rep]] <- est$Phat
    if (!external_reference) {
      Pbar <- Reduce(`+`, aligned_Ps)/length(aligned_Ps)
      reference_P = Pbar
    }
    if (rep == 1) {
      colnames(estimates) <- colnames(reference_P)
      write.csv(estimates[rep, , drop = FALSE], file = glue("{filename}.csv"), row.names = FALSE)
    } else {
      write.table(estimates[rep, , drop = FALSE], file = glue("{filename}.csv"),
                  sep = ",", append = TRUE, col.names = FALSE, row.names = FALSE)
    }
  }
  saveRDS(aligned_Ps, file = glue("{filename}_aligned_Ps.rds"))

  # bootstrap samples have different individuals/orders, can't average Chats
  # instead, back-calculate from bootstrap Phat after the fact
  bootstrap_Phat = Reduce(`+`, aligned_Ps)/length(aligned_Ps)
  bootstrap_Chat = nnlm_wrapper(M, fixed_P = bootstrap_Phat, likelihood = nnlm_likelihood)

  return(list(
    bootstrap_mean = apply(estimates, 2, mean, na.rm = TRUE),
    bootstrap_se = apply(estimates, 2, sd, na.rm = TRUE),
    bootstrap_Phat = bootstrap_Phat,
    bootstrap_Chat = bootstrap_Chat
  ))
}

# ---------

combine_mat <- function(M0, M1, Tr) {
  G <- ncol(M1)
  do.call(cbind, lapply(1:G, function(g) {
    if (Tr[g] == 1) {
      return(M1[,g])
    } else {
      return(M0[,g])
    }
  }))
}

split_dat <- function(M, Tr, prop, FPT = FALSE, PT = 0.5, seed = NULL, force_second = c()) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  G = ncol(M)
  n = round(prop*G) # number in first part
  idx <- 1:G

  if (FPT) {
    # fix proportion treated in part 1
    n1 = round(PT * n)
    n0 = n - n1

    avail_1 <- which(Tr == 1 & !(idx %in% force_second))
    avail_0 <- which(Tr == 0 & !(idx %in% force_second))  

    if (n1 > length(avail_1)) {
      n1 = length(avail_1) - 1
      n0 = round(n1/PT-n1)
      warning(glue("Had to reduce prop to {(n1+n0)/G} to accommodate PT= {PT}"))
    }

    if (n0 > length(avail_0)) {
      n0 = length(avail_0) - 1
      n1 = round(n0/(1-PT)-n0)
      warning(glue("Had to reduce prop to {(n1+n0)/G} to accommodate PT= {PT}"))
      if ((n1 > length(avail_1))) {
        n1 = length(avail_1) - 1
        warning(glue("Can't accommodate PT = {PT}. Using prop {(n1+n0)/G} and PT {n1/(n1+n0)}."))
      
      }
    }

    idx_first = c(
      sample(avail_1, n1),
      sample(avail_0, n0)
    )
  } else {
    # random split
    idx_first = sample(which(!(idx %in% force_second)), size = n)
  }

  return(idx_first)
}

expand <- function(nmf_res_sub, problem_rows) {
  What_sub <- nmf_res_sub@fit@W

  # put 0 back in for problem rows
  What <- matrix(0, nrow = nrow(What_sub) + sum(problem_rows), ncol = ncol(What_sub))
  What[!problem_rows,] <- What_sub

  stopifnot(all((rowSums(What) == 0) == problem_rows))
  nmf_res_sub@fit@W = What
  return(nmf_res_sub)
}

poisson_nnlm <- function(M, fixed_P, maxiter = 1000, tol = 1e-5) {
  K <- nrow(M); G <- ncol(M); N <- ncol(fixed_P)

  # deterministic initialization of equal distribution of mutations across signatures
  C <- matrix(colSums_wrapper(M)/N, nrow = N, ncol = G, byrow = TRUE)

  eps <- 1e-10
  C_prev <- C
  for (iter in 1:maxiter) {
    numerator <- t(fixed_P) %*% (M / fixed_P%*%C)
    denominator <- matrix(colSums_wrapper(fixed_P), nrow = N, ncol = G)
    C <- C * numerator / (denominator + eps)
    delta <- sum(abs(C - C_prev)) / (sum(abs(C_prev)) + eps)
    if (delta < tol) {
      break
    }

    C_prev <- C
  }
  
  return(C)
}

nnlm_wrapper <- function(M, fixed_P, likelihood = 'poisson', bayesian = FALSE, filename = NULL) {
  if (likelihood == 'gaussian') {
    if (bayesian) {
      nnlm_res <- bayesNMF(M, rank = ncol(fixed_P), likelihood = 'normal', fixed = list(P = fixed_P), file = filename)
      C <- nnlm_res$MAP$E
    } else {
      C <- do.call(cbind, lapply(1:ncol(M), function(g) {
        nnls::nnls(A = fixed_P, b = M[,g])$x
      }))
    }
  } else if (likelihood == 'poisson') {
    if (bayesian) {
      nnlm_res <- bayesNMF(M, rank = ncol(fixed_P), fixed = list(P = fixed_P), file = filename)
      C <- nnlm_res$MAP$E
    } else {
      C <- poisson_nnlm(M, fixed_P)
    }
  } else {
    stop(glue("likelihood {likelihood} not implemented for nnlm_wrapper"))
  }
  return(C)
}

nmf_wrapper <- function(M, rank, nrun, seed, method = 'brunet', reference_P = NULL, bayesian = FALSE, filename = NULL) {
  if (bayesian) {
    nmf_res <- bayesNMF(M, rank = rank, file = filename)
    nmf_res <- extract_nmf_info(nmf_res, reference_P, bayesian = bayesian)
  } else{ 
    # problem when one row is all 0
    problem_rows <- rowSums(M) == 0
    if (sum(problem_rows) == 0){
      nmf_res <- NMF::nmf(M, rank = rank, nrun = nrun, seed = seed, method = method)
    } else {
      M_sub <- M[!problem_rows,]
      nmf_res_sub <- NMF::nmf(M_sub, rank = rank, nrun = nrun, seed = seed, method = method)
      nmf_res <- expand(nmf_res_sub, problem_rows)
    }
  }
  return(nmf_res)
}

extract_nmf_info <- function(nmf_res, reference_P = NULL, bayesian = FALSE) {
  if (bayesian) {
    Phat <- nmf_res$MAP$P
    Chat <- nmf_res$MAP$E
  } else {
    Phat <- nmf_res@fit@W
    Chat <- nmf_res@fit@H
  }
  
  N <- ncol(Phat)

  # rescale so E is on the scale of mutation counts
  Chat <- sweep(Chat, 1, colSums(Phat), `*`)
  Phat <- sweep(Phat, 2, colSums(Phat), `/`)

  # reorder estimated signatures to match reference
  minsim = NA
  reassigned = NA
  if (!is.null(reference_P)) {
    sim <- bayesNMF::pairwise_sim(reference_P, Phat, name2 = 'est')
    colnames(Phat) <- paste0('est', 1:N)
    rownames(Chat) <- paste0('est', 1:N)
    reassigned <- bayesNMF::assign_signatures(sim)
    minsim = min(diag(reassigned))
    Phat <- Phat[, colnames(reassigned)]
    Chat <- Chat[colnames(reassigned), ]

    colnames(Phat) <- rownames(reassigned)
    rownames(Chat) <- rownames(reassigned)
  } else {
    print("no reference")
  }

  return(list(
    P = Phat,
    C = Chat,
    minsim = minsim,
    reassigned = reassigned
  ))
}

colSums_wrapper <- function(mat) {
  # allows a vector to be passed in (returns itself)
  # otherwise row means
  if (!('matrix' %in% class(mat))) {
    return(mat)
  }
  return(colSums(mat))
}

rowMeans_wrapper <- function(mat) {
  # allows a vector to be passed in (returns itself)
  # otherwise row means
  if (!('matrix' %in% class(mat))) {
    return(mat)
  }
  return(rowMeans(mat))
}

rowMedians_wrapper <- function(mat) {
  # allows a vector to be passed in (returns itself)
  # otherwise row means
  if (!('matrix' %in% class(mat))) {
    return(mat)
  }
  return(apply(mat, 1, median))
}

decompose <- function(psi_M, P) {
  psi_C <- solve(t(P) %*% P) %*% t(P) %*% psi_M
  return(psi_C[,1])
}
