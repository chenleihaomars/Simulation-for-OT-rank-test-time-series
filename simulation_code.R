# -----------------------------------------------------------------------------
# Simulation design
# -----------------------------------------------------------------------------

CASE_INFO <- data.frame(
  case = seq_len(16L),
  label = c(
    "VAR(1), iid Gaussian",
    "VAR(1), iid t5",
    "VAR(1), iid smoothed shifted exponential",
    "Dense VARMA(1,1), correctly fitted",
    "CCC-GARCH(1,1), correctly fitted",
    "VAR(1) with one unit root, iid Gaussian",
    "Dense linear cross-lag dependence",
    "Dense nonlinear conditional mean",
    "Dense conditional scale",
    "Dense nonlinear lag-3 dependence",
    "Noncommuting cross-MA VARMA(1,1), fitted as VAR(1)",
    "Dense smoothed shifted-exponential lag-1 dependence",
    "CCC-GARCH(1,1), nonlinear dependent innovations",
    "VAR(1) with one unit root, conditional-scale innovations",
    "VAR(1), iid shifted exponential (robustness)",
    "Dense shifted-exponential lag-1 dependence (robustness)"
  ),
  role = c(
    rep("size_main", 6L),
    rep("power_main", 8L),
    "size_robustness",
    "power_robustness"
  ),
  # The nonsequential IDs preserve the random streams of earlier code versions.
  stream_id = c(1:4, 13L, 14L, 5:9, 11L, 15L, 16L, 10L, 12L),
  stringsAsFactors = FALSE
)

SIMULATION_SEED <- 205759L

.model_cache <- new.env(parent = emptyenv())
.worker_cfg <- NULL
.worker_refs <- NULL

split_size <- function(n) {
  r <- as.integer(floor(n/(1.2*log(n))))
  c(s = as.integer(n - r), r = r)
}

truncation_level <- function(r) {
  max(1L, as.integer(floor(r^(1 / 6))))
}

set_stream_seed <- function(case, p, n, irep, stream = 0L) {
  z <- (SIMULATION_SEED + 1000003 * case + 10007 * p + 1009 * n +
          9176 * irep + 101 * stream) %% 2147483647
  set.seed(as.integer(z))
}

# -----------------------------------------------------------------------------
# Model matrices and data generation
# -----------------------------------------------------------------------------

dct_matrix <- function(p) {
  Q <- matrix(0, p, p)
  Q[1L, ] <- 1 / sqrt(p)
  if (p > 1L) {
    j <- 2:p
    for (a in j)
      Q[a, ] <- sqrt(2 / p) *
        cos(pi * (a - 1) * (2 * seq_len(p) - 1) / (2 * p))
  }
  t(Q)
}

model_spec <- function(p) {
  key <- as.character(p)
  if (exists(key, envir = .model_cache, inherits = FALSE))
    return(get(key, envir = .model_cache, inherits = FALSE))
  
  Q <- dct_matrix(p)
  var_eig <- seq(0.15, 0.65, length.out = p)
  va_eig <- seq(0.35, 0.80, length.out = p)
  vb_eig <- seq(-0.55, 0.45, length.out = p)
  # One common stochastic trend; cointegration rank p - 1.
  coint_eig <- c(1, seq(0.30, 0.65, length.out = p - 1L))
  # Constant-correlation matrix for the CCC-GARCH designs.
  R <- matrix(0.35, p, p)
  diag(R) <- 1
  
  spec <- list(
    Q = Q,
    var_eig = var_eig,
    coint_eig = coint_eig,
    va_eig = va_eig,
    vb_eig = vb_eig,
    garch_par = c(omega = 0.05, alpha = 0.08, beta = 0.87),
    garch_chol = chol(R),
    true_a = 0.65,
    true_b = 0.60
  )
  assign(key, spec, envir = .model_cache)
  spec
}

recursive_filter <- function(U, ar) {
  Y <- U
  for (j in seq_len(ncol(U)))
    Y[, j] <- as.numeric(stats::filter(
      U[, j], filter = ar[j], method = "recursive", init = 0
    ))
  Y
}

innov_iid <- function(N, p, type, smooth_exp_tau = 0.25) {
  switch(
    type,
    gaussian = matrix(rnorm(N * p), N, p),
    t5 = {
      df <- 5
      Z <- matrix(rnorm(N * p), N, p)
      Z * sqrt((df - 2) / rchisq(N, df))
    },
    exp = matrix(rexp(N * p) - 1, N, p),
    smooth_exp = {
      E <- matrix(rexp(N * p) - 1, N, p)
      Z <- matrix(rnorm(N * p), N, p)
      (E + smooth_exp_tau * Z) / sqrt(1 + smooth_exp_tau^2)
    },
    stop("Unknown innovation type: ", type)
  )
}

innov_skew_lag1 <- function(N, p, type, smooth_exp_tau = 0.25) {
  # Dense MA(1)-type dependence with identity marginal covariance.
  # The smoothed version is the regular main design; the unsmoothed
  # shifted-exponential version is retained as a robustness design.
  B <- model_spec(p)$Q
  eta <- innov_iid(N + 1L, p, type, smooth_exp_tau)
  current <- eta[2:(N + 1L), , drop = FALSE]
  lag1 <- eta[1:N, , drop = FALSE]
  (current + lag1 %*% t(B)) / sqrt(2)
}

innov_dependent <- function(N, p, type) {
  sp <- model_spec(p)
  B <- sp$Q
  C <- sp$Q[, c(2:p, 1L), drop = FALSE]
  E <- matrix(rnorm((N + 4L) * p), N + 4L, p)
  current <- E[5:(N + 4L), , drop = FALSE]
  lag1 <- E[4:(N + 3L), , drop = FALSE]
  lag3 <- E[2:(N + 1L), , drop = FALSE]
  g <- function(x) (x^2 - 1) / sqrt(2)

  switch(
    type,
    linear = {
      theta <- 1
      (current + theta * (lag1 %*% t(B))) / sqrt(1 + theta^2)
    },
    nonlinear = {
      theta <- 0.9
      (current + theta * (g(lag1 %*% t(C)) %*% t(B))) /
        sqrt(1 + theta^2)
    },
    scale = {
      rho <- 0.8
      d <- exp(0.5 * rho * (lag1 %*% t(C)) - 0.25 * rho^2)
      (d * current) %*% t(B)
    },
    lag3 = {
      theta <- 1
      (current + theta * (g(lag3 %*% t(C)) %*% t(B))) /
        sqrt(1 + theta^2)
    },
    stop("Unknown dependent innovation type: ", type)
  )
}

simulate_var <- function(n, eps) {
  # Latent diagonal filtering is algebraically equivalent to the dense
  # VAR recursion with coefficient Q diag(var_eig) Q'.
  sp <- model_spec(ncol(eps))
  E <- eps %*% sp$Q
  Y <- recursive_filter(E, sp$var_eig)
  X <- Y %*% t(sp$Q)
  N <- nrow(X)
  X[(N - n + 1L):N, , drop = FALSE]
}

simulate_varma <- function(n, p, burnin) {
  sp <- model_spec(p)
  N <- n + burnin
  E <- matrix(rnorm(N * p), N, p) %*% sp$Q
  U <- E
  if (N > 1L)
    U[2:N, ] <- U[2:N, , drop = FALSE] +
    sweep(E[1:(N - 1L), , drop = FALSE], 2,
          sp$true_b * sp$vb_eig, "*")
  Y <- recursive_filter(U, sp$true_a * sp$va_eig)
  X <- Y %*% t(sp$Q)
  X[(burnin + 1L):N, , drop = FALSE]
}

# Unit-root designs start from X_0 = 0 and use no burn-in.
simulate_cointegrated_var <- function(eps) {
  sp <- model_spec(ncol(eps))
  recursive_filter(eps %*% sp$Q, sp$coint_eig) %*% t(sp$Q)
}

simulate_garch <- function(n, eps) {
  p <- ncol(eps)
  sp <- model_spec(p)
  eta <- eps %*% sp$garch_chol
  X <- matrix(0, nrow(eta), p)
  h <- rep(1, p)
  garch <- sp$garch_par
  for (t in seq_len(nrow(X))) {
    if (t > 1L)
      h <- garch["omega"] + garch["alpha"] * X[t - 1L, ]^2 +
        garch["beta"] * h
    X[t, ] <- sqrt(h) * eta[t, ]
  }
  tail(X, n)
}

# Signed cyclic permutation used in case 11. In row-vector notation,
# z %*% R = (-z_p, z_1, ..., z_{p-1}); this is R_p z in the paper's
# column-vector notation. The matrix is orthogonal and
# does not commute with diag(sp$va_eig) when the AR eigenvalues differ.
signed_cyclic_shift <- function(p) {
  R <- matrix(0, p, p)
  if (p > 1L)
    R[cbind(seq_len(p - 1L), 2:p)] <- 1
  R[p, 1L] <- -1
  R
}

# Noncommuting cross-MA VARMA(1,1) alternative:
#   Y_t = D Y_{t-1} + e_t + theta e_{t-1} R,
# followed by the dense rotation X_t = Y_t Q'. Every MA singular value
# equals theta, so the omitted-MA signal remains present in every direction.
simulate_cross_varma <- function(n, p, burnin, theta = 0.70) {
  sp <- model_spec(p)
  N <- n + burnin
  E <- matrix(rnorm(N * p), N, p)
  U <- E
  
  if (N > 1L) {
    R <- signed_cyclic_shift(p)
    U[2:N, ] <- U[2:N, , drop = FALSE] +
      theta * (E[1:(N - 1L), , drop = FALSE] %*% R)
  }
  
  Y <- recursive_filter(U, sp$true_a * sp$va_eig)
  X <- Y %*% t(sp$Q)
  X[(burnin + 1L):N, , drop = FALSE]
}

generate_case <- function(case, n, p, burnin, smooth_exp_tau = 0.25) {
  N <- n + burnin
  result <- function(X, fit = "var") list(X = X, fit = fit)

  switch(
    as.character(case),
    "1" = result(simulate_var(n, innov_iid(N, p, "gaussian"))),
    "2" = result(simulate_var(n, innov_iid(N, p, "t5"))),
    "3" = result(simulate_var(
      n, innov_iid(N, p, "smooth_exp", smooth_exp_tau)
    )),
    "4" = result(simulate_varma(n, p, burnin), "varma"),
    "5" = result(simulate_garch(n, innov_iid(N, p, "gaussian")), "garch"),
    "6" = result(simulate_cointegrated_var(innov_iid(n, p, "gaussian"))),
    "7" = result(simulate_var(n, innov_dependent(N, p, "linear"))),
    "8" = result(simulate_var(n, innov_dependent(N, p, "nonlinear"))),
    "9" = result(simulate_var(n, innov_dependent(N, p, "scale"))),
    "10" = result(simulate_var(n, innov_dependent(N, p, "lag3"))),
    "11" = result(simulate_cross_varma(n, p, burnin, theta = 0.70)),
    "12" = result(simulate_var(
      n, innov_skew_lag1(N, p, "smooth_exp", smooth_exp_tau)
    )),
    "13" = result(
      simulate_garch(n, innov_dependent(N, p, "nonlinear")), "garch"
    ),
    "14" = result(simulate_cointegrated_var(
      innov_dependent(n, p, "scale")
    )),
    "15" = result(simulate_var(
      n, innov_iid(N, p, "exp", smooth_exp_tau)
    )),
    "16" = result(simulate_var(
      n, innov_skew_lag1(N, p, "exp", smooth_exp_tau)
    )),
    stop("Unknown simulation case: ", case)
  )
}

# -----------------------------------------------------------------------------
# Training-block estimation and test-block residuals
# -----------------------------------------------------------------------------

fit_var_test_residuals <- function(X, s) {
  X0 <- X[1:(s - 1L), , drop = FALSE]
  X1 <- X[2:s, , drop = FALSE]
  # The unpenalized Gaussian QMLE equals multivariate least squares.
  # Pivoted QR avoids both the normal equations and ridge regularization.
  B <- qr.solve(X0, X1, tol = 1e-10)
  X[(s + 1L):nrow(X), , drop = FALSE] -
    X[s:(nrow(X) - 1L), , drop = FALSE] %*% B
}

varma_latent_residuals <- function(Y, a, b, sp) {
  n <- nrow(Y)
  U <- matrix(0, n, ncol(Y))
  if (n > 1L)
    U[2:n, ] <- Y[2:n, , drop = FALSE] -
    sweep(Y[1:(n - 1L), , drop = FALSE], 2, a * sp$va_eig, "*")
  recursive_filter(U, -b * sp$vb_eig)
}

fit_varma_test_residuals <- function(X, s) {
  sp <- model_spec(ncol(X))
  Y <- X %*% sp$Q
  Ytr <- Y[1:s, , drop = FALSE]
  
  objective <- function(par) {
    E <- varma_latent_residuals(Ytr, par[1L], par[2L], sp)
    mean(rowSums(E[-1L, , drop = FALSE]^2))
  }
  
  opt <- optim(
    c(0.5, 0.3), objective, method = "L-BFGS-B",
    lower = c(-0.95, -0.95), upper = c(0.95, 0.95),
    control = list(maxit = 60L, factr = 1e8)
  )
  if (opt$convergence != 0L) {
    alt <- nlminb(
      c(0.5, 0.3), objective,
      lower = c(-0.95, -0.95), upper = c(0.95, 0.95),
      control = list(iter.max = 80L, eval.max = 160L)
    )
    if (is.finite(alt$objective) && alt$objective < opt$value)
      opt$par <- alt$par
  }
  
  E <- varma_latent_residuals(Y, opt$par[1L], opt$par[2L], sp)
  E[(s + 1L):nrow(E), , drop = FALSE] %*% t(sp$Q)
}

garch_ab <- function(x) {
  z <- exp(c(0, x) - max(c(0, x)))
  setNames(0.995 * z[2:3] / sum(z), c("alpha", "beta"))
}

garch_variance <- function(X, ab, target) {
  U <- matrix((1 - sum(ab)) * target,
              nrow(X), ncol(X), byrow = TRUE)
  U[1L, ] <- target
  if (nrow(X) > 1L)
    U[2:nrow(X), ] <- U[2:nrow(X), , drop = FALSE] +
      ab["alpha"] * X[1:(nrow(X) - 1L), , drop = FALSE]^2
  recursive_filter(U, rep(ab["beta"], ncol(X)))
}

fit_garch_test_residuals <- function(X, s) {
  # Common persistence parameters with coordinatewise variance targeting.
  Xtr <- X[1:s, , drop = FALSE]
  target <- pmax(colMeans(Xtr^2), 1e-8)
  #The true values are used only as the optimizer’s starting point. 
  #It does not fix the estimator at the truth: optim() is free to move away from these values.
  start <- log(c(0.08, 0.87) / (0.995 - 0.95))
  objective <- function(x) {
    H <- garch_variance(Xtr, garch_ab(x), target)
    mean(log(H[-1L, , drop = FALSE]) +
           Xtr[-1L, , drop = FALSE]^2 / H[-1L, , drop = FALSE])
  }
  opt <- optim(start, objective, method = "BFGS",
               control = list(maxit = 60L, reltol = 1e-7))
  if (opt$convergence != 0L || !is.finite(opt$value)) {
    alt <- nlminb(start, objective,
                  control = list(iter.max = 80L, eval.max = 160L))
    if (is.finite(alt$objective) &&
        (!is.finite(opt$value) || alt$objective < opt$value))
      opt$par <- alt$par
  }
  if (any(!is.finite(opt$par)) || !is.finite(objective(opt$par)))
    stop("GARCH QMLE failed.")
  H <- garch_variance(X, garch_ab(opt$par), target)
  X[(s + 1L):nrow(X), , drop = FALSE] /
    sqrt(H[(s + 1L):nrow(H), , drop = FALSE])
}

fit_test_residuals <- function(obj, s) {
  if (obj$fit == "varma") return(fit_varma_test_residuals(obj$X, s))
  if (obj$fit == "garch") return(fit_garch_test_residuals(obj$X, s))
  fit_var_test_residuals(obj$X, s)
}

# -----------------------------------------------------------------------------
# Tensor-cosine features, OT assignment, and test calibration
# -----------------------------------------------------------------------------

multi_index_spec <- function(p, K) {
  index <- as.matrix(expand.grid(rep(list(0:K), p),
                                 KEEP.OUT.ATTRS = FALSE))
  degree <- rowSums(index)
  index <- index[degree > 0 & degree <= K, , drop = FALSE]
  degree <- rowSums(index)
  omega <- function(q) ifelse(q == 0, 1, 6 / (pi^2 * q^2))
  lambda <- apply(index, 1L, function(z) prod(omega(z))) / (2^p - 1)
  ord <- order(degree, -lambda)
  list(index = index[ord, , drop = FALSE],
       lambda = lambda[ord], d = length(lambda))
}

basis_on_grid <- function(U, spec) {
  r <- nrow(U)
  p <- ncol(U)
  K <- max(spec$index)
  marginal <- vector("list", p)
  for (j in seq_len(p)) {
    A <- matrix(1, r, K + 1L)
    for (q in seq_len(K))
      A[, q + 1L] <- sqrt(2) * cos(pi * q * (U[, j] + 1) / 2)
    marginal[[j]] <- A
  }
  
  B <- matrix(1, r, spec$d)
  for (j in seq_len(p))
    B <- B * marginal[[j]][, spec$index[, j] + 1L, drop = FALSE]
  sweep(B, 2, sqrt(spec$lambda), "*")
}

solve_assignment_score <- function(score, solver) {
  if (solver == "iGraphMatch") {
    m <- iGraphMatch::do_lap(score, method = "lapjv")
  } else if (solver == "TreeDist") {
    m <- TreeDist::LAPJV(max(score) - score)$matching
  } else if (solver == "clue") {
    score <- score - min(score)
    m <- clue::solve_LSAP(score, maximum = TRUE)
  } else {
    stop("Unknown assignment solver.")
  }
  m <- as.integer(m)
  n <- nrow(score)
  if (length(m) != n || anyNA(m) ||
      !identical(sort.int(m), seq_len(n)))
    stop("The assignment solver returned an invalid permutation.")
  m
}

resolve_solver <- function(requested = "auto", allow_slow = FALSE) {
  solver_names <- c("iGraphMatch", "TreeDist", "clue")
  available <- setNames(
    vapply(solver_names, requireNamespace, logical(1L), quietly = TRUE),
    solver_names
  )
  
  if (requested != "auto") {
    if (!requested %in% names(available) || !available[requested])
      stop("Requested solver is unavailable: ", requested)
    return(list(name = requested, benchmark = NA_real_))
  }
  
  fast <- names(available)[available & names(available) != "clue"]
  if (!length(fast)) {
    if (!allow_slow)
      stop(
        "No fast Jonker--Volgenant solver is available. ",
        "Install the binary package 'iGraphMatch' and restart R."
      )
    return(list(name = "clue", benchmark = NA_real_))
  }
  
  set.seed(SIMULATION_SEED)
  Z <- matrix(rnorm(300L * 3L), 300L, 3L)
  G <- matrix(runif(300L * 3L, -1, 1), 300L, 3L)
  score <- tcrossprod(Z, G)
  
  reference <- as.integer(clue::solve_LSAP(
    score - min(score), maximum = TRUE
  ))
  reference_value <- sum(score[cbind(seq_len(nrow(score)), reference)])
  
  elapsed <- setNames(rep(Inf, length(fast)), fast)
  for (z in fast) {
    ans <- tryCatch({
      m <- solve_assignment_score(score, z)
      value <- sum(score[cbind(seq_len(nrow(score)), m)])
      if (!isTRUE(all.equal(value, reference_value, tolerance = 1e-7)))
        stop("objective mismatch")
      median(replicate(
        2L,
        unname(system.time(solve_assignment_score(score, z))["elapsed"])
      ))
    }, error = function(e) Inf)
    elapsed[z] <- ans
  }
  
  if (all(!is.finite(elapsed))) {
    if (!allow_slow)
      stop("Installed fast solvers failed their correctness check.")
    return(list(name = "clue", benchmark = NA_real_))
  }
  name <- names(which.min(elapsed))
  list(name = name, benchmark = unname(elapsed[name]))
}

ot_assignment <- function(Z, G, solver) {
  if (any(!is.finite(Z))) stop("Non-finite fitted residual.")
  # The row and column norm terms are constant over permutations, so
  # squared-Euclidean OT is equivalent to maximizing Z %*% t(G).
  solve_assignment_score(tcrossprod(Z, G), solver)
}

lag_statistics <- function(F, Mvals) {
  r <- nrow(F)
  norm2 <- rowSums(F^2)
  Tmax <- max(Mvals)
  T <- numeric(Tmax)
  
  for (m in seq_len(Tmax)) {
    N <- r - m
    A <- F[seq_len(N), , drop = FALSE]
    B <- F[(m + 1L):(m + N), , drop = FALSE]
    C <- crossprod(A, B)
    T[m] <- (sum(C^2) -
               sum(norm2[seq_len(N)] * norm2[(m + 1L):(m + N)])) / N
  }
  cumsum(T^2)[Mvals]
}

critical_values <- function(lambda, Mvals, ndraw, alpha) {
  w <- as.vector(outer(lambda, lambda))
  tab <- table(signif(w, 14))
  value <- as.numeric(names(tab))
  df <- as.integer(tab)
  Tmax <- max(Mvals)
  V <- matrix(0, ndraw, Tmax)
  
  for (j in seq_along(value))
    V <- V + value[j] *
    matrix(rchisq(ndraw * Tmax, df[j]) - df[j], ndraw, Tmax)
  
  S <- V^2
  if (Tmax > 1L)
    for (m in 2:Tmax) S[, m] <- S[, m] + S[, m - 1L]
  
  vapply(
    Mvals,
    function(m) unname(quantile(S[, m], 1 - alpha, type = 8)),
    numeric(1L)
  )
}

prepare_reference <- function(cfg) {
  refs <- list()
  cache <- list()
  
  for (p in cfg$p) for (n in cfg$n) {
    split <- split_size(n)
    K <- truncation_level(split["r"])
    key <- paste(p, K, sep = "_")
    
    if (is.null(cache[[key]])) {
      spec <- multi_index_spec(p, K)
      set.seed(SIMULATION_SEED + 31L * p + 7919L * K)
      cache[[key]] <- list(
        spec = spec,
        critical = critical_values(
          spec$lambda, cfg$M, cfg$null_draws, cfg$alpha
        )
      )
    }
    
    refs[[paste(p, n, sep = "_")]] <- c(
      list(s = unname(split["s"]), r = unname(split["r"]), K = K),
      cache[[key]]
    )
  }
  refs
}

# -----------------------------------------------------------------------------
# Monte Carlo engine and output
# -----------------------------------------------------------------------------

simulate_replication <- function(p, n, irep, cfg, refs) {
  ref <- refs[[paste(p, n, sep = "_")]]
  r <- ref$r
  s <- ref$s
  nc <- length(cfg$cases)
  nm <- length(cfg$M)
  
  set_stream_seed(0L, p, n, irep, 1L)
  grid <- matrix(runif(r * p, -1, 1), r, p)
  grid_feature <- basis_on_grid(grid, ref$spec)
  
  statistic <- matrix(NA_real_, nc, nm)
  reject <- matrix(NA, nc, nm)
  error <- character(nc)
  
  for (ii in seq_along(cfg$cases)) {
    case <- cfg$cases[ii]
    set_stream_seed(CASE_INFO$stream_id[case], p, n, irep, 2L)
    
    tryCatch({
      tau <- if (is.null(cfg$smooth_exp_tau)) 0.25 else cfg$smooth_exp_tau
      obj <- generate_case(case, n, p, cfg$burnin, tau)
      Z <- fit_test_residuals(obj, s)
      matching <- ot_assignment(Z, grid, cfg$solver)
      statistic[ii, ] <- lag_statistics(
        grid_feature[matching, , drop = FALSE], cfg$M
      )
      reject[ii, ] <- statistic[ii, ] > ref$critical
    }, error = function(e) {
      error[ii] <<- conditionMessage(e)
    })
  }
  
  data.frame(
    case = rep(cfg$cases, each = nm),
    p = p,
    n = n,
    s = s,
    r = r,
    K = ref$K,
    rep = irep,
    M = rep(cfg$M, times = nc),
    statistic = as.vector(t(statistic)),
    critical = rep(ref$critical, times = nc),
    reject = as.vector(t(reject)),
    error = rep(error, each = nm),
    row.names = NULL
  )
}

worker_replication <- function(irep, p, n) {
  simulate_replication(p, n, irep, .worker_cfg, .worker_refs)
}

valid_partial_block <- function(x, p, n, cfg) {
  required <- c("case", "p", "n", "rep", "M", "reject", "error")
  if (!is.data.frame(x) || !all(required %in% names(x)) ||
      any(x$p != p) || any(x$n != n)) return(FALSE)
  
  expected_per_rep <- length(cfg$cases) * length(cfg$M)
  tab <- table(x$rep)
  identical(sort(unique(x$case)), sort(as.integer(cfg$cases))) &&
    identical(sort(unique(x$M)), sort(as.integer(cfg$M))) &&
    all(tab == expected_per_rep) &&
    !anyDuplicated(paste(x$rep, x$case, x$M, sep = ":")) &&
    all(as.integer(names(tab)) %in% seq_len(cfg$reps))
}

format_seconds <- function(x) {
  if (!is.finite(x)) return("?")
  if (x < 60) return(sprintf("%.0f sec", x))
  if (x < 3600) return(sprintf("%.1f min", x / 60))
  sprintf("%.2f hr", x / 3600)
}

run_one_block <- function(p, n, cfg, refs, cluster, block_file) {
  existing <- NULL
  if (cfg$resume && file.exists(block_file)) {
    z <- tryCatch(readRDS(block_file), error = function(e) NULL)
    if (!is.null(z) && valid_partial_block(z, p, n, cfg)) existing <- z
  }
  
  done <- if (is.null(existing)) integer() else sort(unique(existing$rep))
  todo <- setdiff(seq_len(cfg$reps), done)
  if (!length(todo)) return(list(data = existing, elapsed = 0))
  
  r <- refs[[paste(p, n, sep = "_")]]$r
  multiplier <- if (r < 250) 8L else if (r < 700) 4L else cfg$wave_per_worker
  wave_size <- max(1L, cfg$ncores * multiplier)
  started <- proc.time()[3L]
  session_done <- 0L
  data <- existing
  
  while (length(todo)) {
    ids <- head(todo, wave_size)
    parts <- if (is.null(cluster)) {
      lapply(ids, function(z) simulate_replication(p, n, z, cfg, refs))
    } else {
      parallel::parLapplyLB(
        cluster, ids, worker_replication, p = p, n = n
      )
    }
    
    add <- do.call(rbind, parts)
    failed <- is.na(add$reject) | !is.finite(add$statistic)
    if (any(failed)) {
      bad <- add[failed, c("case", "p", "n", "rep", "error")]
      file <- file.path(cfg$outdir, sprintf("failures_p%d_n%d.csv", p, n))
      write.csv(bad, file, row.names = FALSE)
      stop(
        "Simulation failed for p=", p, ", n=", n,
        ". First error: ", bad$error[1L], ". See ", file
      )
    }
    
    data <- if (is.null(data)) add else rbind(data, add)
    saveRDS(data, block_file, compress = FALSE)
    session_done <- session_done + length(ids)
    done <- sort(unique(data$rep))
    todo <- setdiff(seq_len(cfg$reps), done)
    
    elapsed <- proc.time()[3L] - started
    speed <- session_done / elapsed
    eta <- if (speed > 0) length(todo) / speed else Inf
    message(
      sprintf(
        "    replications %d/%d | elapsed %s | ETA %s",
        length(done), cfg$reps, format_seconds(elapsed),
        format_seconds(eta)
      )
    )
  }
  
  data <- data[order(data$rep, data$case, data$M), ]
  saveRDS(data, block_file, compress = FALSE)
  list(
    data = data,
    elapsed = proc.time()[3L] - started
  )
}

summarize_results <- function(raw) {
  key <- paste(raw$case, raw$p, raw$n, raw$M, sep = "_")
  out <- lapply(split(raw, key), function(x) {
    rate <- mean(x$reject)
    data.frame(
      case = x$case[1L],
      label = CASE_INFO$label[x$case[1L]],
      role = CASE_INFO$role[x$case[1L]],
      p = x$p[1L],
      n = x$n[1L],
      s = x$s[1L],
      r = x$r[1L],
      K = x$K[1L],
      M = x$M[1L],
      replications = nrow(x),
      rejection_rate = rate,
      mc_se = sqrt(rate * (1 - rate) / nrow(x))
    )
  })
  ans <- do.call(rbind, out)
  rownames(ans) <- NULL
  ans[order(ans$case, ans$p, ans$n, ans$M), ]
}

write_outputs <- function(summary, raw, cfg, refs, timings, solver_info) {
  design <- do.call(rbind, lapply(names(refs), function(key) {
    z <- refs[[key]]
    id <- strsplit(key, "_", fixed = TRUE)[[1L]]
    data.frame(
      p = as.integer(id[1L]),
      n = as.integer(id[2L]),
      s = z$s,
      r = z$r,
      K = z$K,
      basis_dimension = z$spec$d
    )
  }))
  design <- design[order(design$p, design$n), ]
  settings <- data.frame(
    parameter = names(cfg),
    value = vapply(cfg, paste, character(1L), collapse = ", ")
  )
  solver <- data.frame(
    solver = solver_info$name,
    benchmark_seconds_300 = solver_info$benchmark
  )
  cases <- CASE_INFO[, c("case", "label", "role")]
  is_size <- grepl("^size_", summary$role)
  is_power <- grepl("^power_", summary$role)
  is_robustness <- grepl("_robustness$", summary$role)
  
  writexl::write_xlsx(
    list(
      Summary = summary,
      Size_Main = summary[summary$role == "size_main", ],
      Power_Main = summary[summary$role == "power_main", ],
      Robustness = summary[is_robustness, ],
      Size = summary[is_size, ],
      Power = summary[is_power, ],
      Size_Robustness = summary[summary$role == "size_robustness", ],
      Power_Robustness = summary[summary$role == "power_robustness", ],
      Design = design,
      Cases = cases,
      Timings = timings,
      Solver = solver,
      Settings = settings
    ),
    file.path(cfg$outdir, "simulation_results.xlsx")
  )
  write.csv(summary, file.path(cfg$outdir, "summary_results.csv"),
            row.names = FALSE)
  saveRDS(raw, file.path(cfg$outdir, "raw_results.rds"),
          compress = "gzip")
}

run_experiment <- function(cfg, core_file) {
  cfg$seed <- SIMULATION_SEED
  
  if (any(!cfg$cases %in% CASE_INFO$case))
    stop("cfg$cases contains an unknown case number.")
  if (!is.null(cfg$smooth_exp_tau) &&
      (!is.finite(cfg$smooth_exp_tau) || cfg$smooth_exp_tau <= 0))
    stop("cfg$smooth_exp_tau must be a positive finite number.")
  
  if (!dir.exists(cfg$outdir) &&
      !dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE))
    stop("Cannot create output directory: ", cfg$outdir)
  probe <- file.path(cfg$outdir, ".write_test")
  if (!file.create(probe)) stop("The code folder is not writable.")
  unlink(probe)
  
  solver_info <- resolve_solver(cfg$solver, cfg$allow_slow_solver)
  cfg$solver <- solver_info$name
  message(
    "Assignment solver: ", cfg$solver,
    if (is.finite(solver_info$benchmark))
      paste0(" (300 x 300 benchmark: ",
             sprintf("%.3f", solver_info$benchmark), " sec)") else ""
  )
  
  message("Running serial preflight for all configured cases...")
  pre <- cfg
  pre$null_draws <- 2000L
  pre_refs <- prepare_reference(pre)
  chk <- simulate_replication(cfg$p[1L], cfg$n[1L], 1L, pre, pre_refs)
  if (anyNA(chk$reject))
    stop("Preflight failed: ", paste(unique(chk$error[nzchar(chk$error)]),
                                     collapse = " | "))
  message("Preflight passed. Preparing Gaussian-chaos critical values...")
  refs <- prepare_reference(cfg)
  
  code_id <- substr(unname(tools::md5sum(core_file)), 1L, 10L)
  cfg_file <- tempfile(fileext = ".rds")
  on.exit(unlink(cfg_file), add = TRUE)
  checkpoint_cfg <- cfg[setdiff(names(cfg), c("ncores", "resume", "outdir"))]
  saveRDS(checkpoint_cfg, cfg_file)
  cfg_id <- substr(unname(tools::md5sum(cfg_file)), 1L, 10L)
  block_dir <- file.path(
    cfg$outdir,
    paste0("blocks_", code_id, "_", cfg_id, "_", cfg$solver,
           "_R", cfg$reps)
  )
  dir.create(block_dir, recursive = TRUE, showWarnings = FALSE)
  
  cluster <- NULL
  if (cfg$ncores > 1L) {
    cluster <- parallel::makeCluster(cfg$ncores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    
    lib_paths <- .libPaths()
    code_path <- normalizePath(core_file, winslash = "/", mustWork = TRUE)
    parallel::clusterExport(
      cluster, c("lib_paths", "code_path"), envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      .libPaths(lib_paths)
      Sys.setenv(
        OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1"
      )
      source(code_path, local = .GlobalEnv)
      NULL
    })
    
    parallel::clusterExport(
      cluster, c("cfg", "refs"), envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      assign(".worker_cfg", cfg, envir = .GlobalEnv)
      assign(".worker_refs", refs, envir = .GlobalEnv)
      rm(cfg, refs)
      NULL
    })
  }
  
  design <- expand.grid(p = cfg$p, n = cfg$n, KEEP.OUT.ATTRS = FALSE)
  blocks <- vector("list", nrow(design))
  timings <- vector("list", nrow(design))
  
  for (i in seq_len(nrow(design))) {
    p <- design$p[i]
    n <- design$n[i]
    message(sprintf("[%02d/%02d] p=%d, n=%d", i, nrow(design), p, n))
    file <- file.path(block_dir, sprintf("block_p%d_n%d.rds", p, n))
    z <- run_one_block(p, n, cfg, refs, cluster, file)
    blocks[[i]] <- z$data
    timings[[i]] <- data.frame(
      p = p, n = n, elapsed_seconds = z$elapsed,
      elapsed = format_seconds(z$elapsed)
    )
  }
  
  raw <- do.call(rbind, blocks)
  rownames(raw) <- NULL
  summary <- summarize_results(raw)
  timing_table <- do.call(rbind, timings)
  write_outputs(summary, raw, cfg, refs, timing_table, solver_info)
  message("Finished. Spreadsheet: ",
          file.path(cfg$outdir, "simulation_results.xlsx"))
  summary
}