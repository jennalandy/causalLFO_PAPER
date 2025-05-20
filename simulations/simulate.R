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

simulate_BA <- function(
    P, C0_dist,
    ATE_C, C0_sd, C1_sd, pi,
    G = 100
) {
  K = nrow(P)
  N = ncol(P)
  idx = sample(1:ncol(C0_dist), G, replace = TRUE)

  # sample with replacement from untreated C sampling distribution
  sampled = C0_dist[,idx]

  # C0 is sampled + Gaussian noise, C1 is sampled + ATE + Gaussian noise
  C0 = do.call(cbind, lapply(1:G, function(g) {
    sampled[,g] + rnorm(N, 0, sd = C0_sd)
  }))
  C1 = do.call(cbind, lapply(1:G, function(g) {
    sampled[,g] + rnorm(N, ATE_C, sd = C1_sd)
  }))

  # ensure non-negativity in both
  C1[C1<0] = 0
  C0[C0<0] = 0

  # generate observed data M from latents C and P
  M0 = do.call(cbind, lapply(1:G, function(g) {
    rpois(K, P%*%C0[,g])
  }))
  M1 = do.call(cbind, lapply(1:G, function(g) {
    rpois(K, P%*%C1[,g])
  }))

  # simulate treatment assignment vector
  Tr = rbinom(G, size = 1, prob = pi)

  dat = list(
    params = list(
      ATE_C = ATE_C,
      C0_sd = C0_sd,
      C1_sd = C1_sd,
      pi = pi,
      P = P,
      C0_dist = C0_dist
    ),
    C0 = C0,
    C1 = C1,
    M0 = M0,
    M1 = M1,
    Tr = Tr,
    M = combine_mat(M0, M1, Tr),
    C = combine_mat(C0, C1, Tr)
  )
}
