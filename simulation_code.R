# OT-rank tests: data generation, fitting, statistics, and parallel execution.
# Matrices store observations in rows and scalar coordinates in columns.

# -----------------------------------------------------------------------------
# Simulation design
# -----------------------------------------------------------------------------

SIMULATION_SEED <- 205759L
# Map F1--F14 to the original stream IDs to reproduce the reported results.
FITTED_STREAM_ID <- c(1:4, 13L, 14L, 5:9, 11L, 15L, 16L)

.model_cache <- new.env(parent = emptyenv())

split_size <- function(n) {
  r <- as.integer(floor(n / (1.2 * log(n))))
  c(s = as.integer(n - r), r = r)
}

truncation_level <- function(r) {
  max(1L, as.integer(floor(r^(1 / 6))))
}

set_stream_seed <- function(case, p, n, irep, stream = 0L) {
  z <- (SIMULATION_SEED + 1000003 * case + 10007 * p + 1009 * n +
          9176 * irep + 101 * stream) %% 2147483647
  set.seed(as.integer(z), kind = "Mersenne-Twister",
           normal.kind = "Inversion", sample.kind = "Rejection")
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

# F1--F14 in the paper; F1--F6 are null designs.
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
    # F11 deliberately omits the moving-average component from the fitted VAR.
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
    stop("Unknown simulation case: ", case)
  )
}

# Directly observed designs O1--O7.
observed_case <- function(case, n, p) {
  if (case == 1L) return(innov_iid(n, p, "gaussian"))
  if (case == 2L) return(innov_iid(n, p, "exp"))
  if (case == 3L) return(innov_dependent(n, p, "lag3"))

  if (case == 4L) {
    # Independent pairs with a dependent sign product in each length-three block.
    phase <- sample.int(3L, 1L) - 1L
    blocks <- ceiling((n + 2L) / 3L)
    a <- matrix(sample(c(-1, 1), blocks * p, TRUE), blocks, p)
    b <- matrix(sample(c(-1, 1), blocks * p, TRUE), blocks, p)
    signs <- matrix(0, 3L * blocks, p)
    signs[seq(1L, 3L * blocks, 3L), ] <- a
    signs[seq(2L, 3L * blocks, 3L), ] <- b
    signs[seq(3L, 3L * blocks, 3L), ] <- a * b
    return(signs[seq_len(n) + phase, , drop = FALSE] *
             matrix(runif(n * p), n, p))
  }

  if (case %in% c(5L, 6L)) {
    u <- runif(n + 1L, -1, 1)
    X <- matrix(runif(n * p, -1, 1), n, p)
    X[, 1L] <- u[-1L]
    if (case == 5L) {
      # Alternating signs cancel the pooled adjacent-pair signal.
      X[, 2L] <- (-1)^seq_len(n) * sign(u[-(n + 1L)]) *
        sample(c(-1, 1), n, TRUE, c(0.1, 0.9)) * runif(n)
    } else {
      # Rejection sampling from f(v | u) = {1 + c(u) sin(pi v / 2)} / 2.
      c <- 0.8 * cos(pi * u[-(n + 1L)])
      pending <- seq_len(n)
      while (length(pending)) {
        v <- runif(length(pending), -1, 1)
        accept <- runif(length(pending)) * (1 + abs(c[pending])) <=
          1 + c[pending] * sin(pi * v / 2)
        X[pending[accept], 2L] <- v[accept]
        pending <- pending[!accept]
      }
    }
    return(X)
  }

  # O7: a persistent latent state switches the within-vector correlation.
  state <- sample(c(-1, 1), 1L) *
    cumprod(c(1, sample(c(-1, 1), n - 1L, TRUE, c(0.1, 0.9))))
  X <- matrix(rnorm(n * p), n, p)
  X[, 2L] <- 0.8 * state * X[, 1L] + 0.6 * X[, 2L]
  X
}

# -----------------------------------------------------------------------------
# Training-block estimation and test-block residuals
# -----------------------------------------------------------------------------

fit_var_residuals <- function(X, s) {
  X0 <- X[1:(s - 1L), , drop = FALSE]
  X1 <- X[2:s, , drop = FALSE]
  # The unpenalized Gaussian QMLE equals multivariate least squares.
  # Pivoted QR avoids both the normal equations and ridge regularization.
  B <- qr.solve(X0, X1, tol = 1e-10)
  n <- nrow(X)
  rbind(rep(0, ncol(X)), X[2:n, , drop = FALSE] -
          X[1:(n - 1L), , drop = FALSE] %*% B)
}

varma_latent_residuals <- function(Y, a, b, sp) {
  n <- nrow(Y)
  U <- matrix(0, n, ncol(Y))
  if (n > 1L)
    U[2:n, ] <- Y[2:n, , drop = FALSE] -
    sweep(Y[1:(n - 1L), , drop = FALSE], 2, a * sp$va_eig, "*")
  recursive_filter(U, -b * sp$vb_eig)
}

fit_varma_residuals <- function(X, s) {
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
  E %*% t(sp$Q)
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

fit_garch_residuals <- function(X, s) {
  # Common persistence parameters with coordinatewise variance targeting.
  Xtr <- X[1:s, , drop = FALSE]
  target <- pmax(colMeans(Xtr^2), 1e-8)
  # The initial values do not constrain the optimizer.
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
  X / sqrt(H)
}

fit_residuals <- function(obj, s) {
  if (obj$fit == "varma") return(fit_varma_residuals(obj$X, s))
  if (obj$fit == "garch") return(fit_garch_residuals(obj$X, s))
  fit_var_residuals(obj$X, s)
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

ot_assignment <- function(Z, G) {
  # For equal-size samples, minimizing squared distances maximizes this score.
  score <- tcrossprod(Z, G)
  as.integer(clue::solve_LSAP(score - min(score), maximum = TRUE))
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

# Deterministic Gaussian-chaos cutoff. The feature columns already include
# sqrt(lambda); the reference is sum_m {sum_j w_j (chi^2_df[j] - df[j])}^2.
critical_values <- function(lambda, Mvals, alpha = 0.05, ngrid = 2^15) {
  tab <- table(signif(as.vector(outer(lambda, lambda)), 14))
  w <- as.numeric(names(tab))
  df <- as.integer(tab)
  scale <- sqrt(2 * sum(df * w^2))
  w <- w / scale

  # Fourier inversion for the centered chaos, scaled to variance one.
  lower <- -sum(df * w) - 1
  upper <- sqrt(60) + 60 * max(w) + 1
  dx <- (upper - lower) / ngrid
  x <- lower + (0:(ngrid - 1)) * dx
  t <- 2 * pi * c(0:(ngrid / 2 - 1), (-ngrid / 2):(-1)) /
    (upper - lower)
  logcf <- -1i * t * sum(df * w)
  for (j in seq_along(w))
    logcf <- logcf - df[j] / 2 * log(1 - 2i * t * w[j])
  density <- pmax(Re(fft(exp(logcf - 1i * t * lower))) /
                    (ngrid * dx), 0)
  cdf <- c(0, cumsum((head(density, -1) + tail(density, -1)) * dx / 2))
  F <- approxfun(x, cdf / tail(cdf, 1), yleft = 0, yright = 1)

  # Bin probabilities for the square; zero padding prevents circular aliasing.
  h <- max(Mvals) / (alpha * ngrid)
  edge <- sqrt((0:ngrid) * h)
  probability <- diff(F(edge) - F(-edge))
  nfft <- 2^ceiling(log2(max(Mvals) * ngrid))
  transform <- fft(c(probability, rep(0, nfft - ngrid)))
  vapply(Mvals, function(m) {
    mass <- pmax(Re(fft(transform^m, inverse = TRUE)) / nfft, 0)
    cumulative <- cumsum(mass)
    j <- which(cumulative >= 1 - alpha)[1]
    # Half-cell correction accounts for summing midpoint-discretized squares.
    quantile <- (j - 1 + (m + 1) / 2 -
                   (cumulative[j] - (1 - alpha)) / mass[j]) * h
    scale^2 * quantile
  }, numeric(1))
}

# Group block coordinates by their interior frequencies and endpoint degree.
# Setting lambda = 1 lets basis_on_grid return the unweighted cosine features.
block_spec <- function(p, L, D) {
  feature <- multi_index_spec(p, D - 1L)
  feature$lambda[] <- 1
  degree <- rowSums(feature$index)
  keys <- c(paste(rep(0, p), collapse = ","),
            apply(feature$index, 1L, paste, collapse = ","))
  interior <- vector("list", L)
  J <- 0L

  for (ell in seq_len(L)) {
    indices <- matrix(0L, 1L, ell - 1L)
    degrees <- 0L
    variables <- p * (ell - 1L)
    if (variables > 0L && D > 2L) {
      for (k in seq_len(D - 2L)) {
        factors <- combn(variables + k - 1L, k) - (seq_len(k) - 1L)
        frequencies <- matrix(0L, ncol(factors), variables)
        for (j in seq_len(k)) {
          pos <- cbind(seq_len(nrow(frequencies)), factors[j, ])
          frequencies[pos] <- frequencies[pos] + 1L
        }
        index <- matrix(0L, nrow(frequencies), ell - 1L)
        for (j in seq_len(ell - 1L)) {
          q <- frequencies[, (j - 1L) * p + seq_len(p), drop = FALSE]
          index[, j] <- match(apply(q, 1L, paste, collapse = ","), keys) - 1L
        }
        indices <- rbind(indices, index)
        degrees <- c(degrees, rep(k, nrow(index)))
      }
    }
    interior[[ell]] <- list(index = indices, degree = degrees)
    for (k in degrees) for (j in seq_len(D - k - 1L))
      J <- J + sum(degree == j) * sum(degree <= D - k - j)
  }
  list(feature = feature, degree = degree, interior = interior,
       L = L, D = D, J = J)
}

# F is basis_on_grid(W, spec$feature), in the observed time order.
block_statistic <- function(F, spec) {
  n <- nrow(F)
  V <- 0
  for (ell in seq_len(spec$L)) {
    t <- seq_len(n - ell)
    left <- lapply(seq_len(spec$D - 1L), function(j)
      F[t, spec$degree == j, drop = FALSE])
    right <- lapply(seq_len(spec$D - 1L), function(j)
      F[t + ell, spec$degree <= j, drop = FALSE])
    interior <- spec$interior[[ell]]
    for (a in seq_len(nrow(interior$index))) {
      h <- rep(1, length(t))
      for (j in which(interior$index[a, ] > 0L))
        h <- h * F[t + j, interior$index[a, j]]
      for (j in seq_len(spec$D - interior$degree[a] - 1L)) {
        C <- crossprod(left[[j]] * h,
                       right[[spec$D - interior$degree[a] - j]])
        V <- max(V, max(C^2) / length(t))
      }
    }
  }
  V
}

# -----------------------------------------------------------------------------
# References and replication-level statistics
# -----------------------------------------------------------------------------

prepare_observed_reference <- function(cfg) {
  refs <- list()
  cache <- list()
  for (p in cfg$p) for (n in cfg$n) {
    K <- truncation_level(n)
    L <- max(3L, floor(log2(n) / 3))
    D <- max(4L, floor(log2(n) / 4))
    key <- paste(p, K, L, D, sep = "_")
    if (is.null(cache[[key]])) {
      pair <- multi_index_spec(p, K)
      block <- block_spec(p, L, D)
      columns <- match(apply(pair$index, 1, paste, collapse = ","),
                       apply(block$feature$index, 1, paste, collapse = ","))
      critical <- c(critical_values(pair$lambda, cfg$M, cfg$alpha),
                    qnorm((1 + (1 - cfg$alpha)^(1 / block$J)) / 2)^2)
      cache[[key]] <- list(K = K, pair = pair, block = block,
                           columns = columns, critical = critical)
    }
    refs[[paste(p, n, sep = "_")]] <- cache[[key]]
  }
  refs
}

observed_statistics <- function(F, ref, Mvals) {
  pair <- sweep(F[, ref$columns, drop = FALSE],
                2, sqrt(ref$pair$lambda), "*")
  c(lag_statistics(pair, Mvals), block_statistic(F, ref$block))
}

simulate_observed_replication <- function(p, n, irep, cfg, refs) {
  ref <- refs[[paste(p, n, sep = "_")]]
  set_stream_seed(100L, p, n, irep, 1L)
  grid <- matrix(runif(n * p, -1, 1), n, p)
  features <- basis_on_grid(grid, ref$block$feature)
  statistic <- t(vapply(cfg$cases, function(case) {
    set_stream_seed(100L + case, p, n, irep, 2L)
    X <- observed_case(case, n, p)
    matching <- ot_assignment(X, grid)
    observed_statistics(features[matching, , drop = FALSE], ref, cfg$M)
  }, numeric(length(cfg$M) + 1L)))
  tests <- c(rep("pairwise", length(cfg$M)), "block")
  data.frame(case = rep(cfg$cases, each = length(tests)), p = p, n = n, rep = irep,
             test = rep(tests, length(cfg$cases)),
             M = rep(c(cfg$M, 0L), length(cfg$cases)),
             statistic = as.vector(t(statistic)),
             critical = rep(ref$critical, length(cfg$cases)),
             reject = as.vector(t(sweep(statistic, 2, ref$critical, ">"))))
}

# The mask retains pairs of nonconstant cosines with total degree at most 3.
prepare_fitted_reference <- function(cfg) {
  refs <- list()
  for (p in cfg$p) for (n in cfg$n) {
    sizes <- split_size(n)
    s <- unname(sizes["s"])
    r <- unname(sizes["r"])
    spec <- multi_index_spec(p, 2L)
    degree <- rowSums(spec$index)
    pair <- outer(degree, degree, "+") <= 3L
    discard <- max(2L, ceiling(s^(1 / 4)))
    chunk_size <- max(32L, floor(s^(3 / 4)))
    J <- max(1L, ceiling((s - discard) / chunk_size))
    chunks <- split(seq.int(discard + 1L, s),
                    ceiling(seq_len(s - discard) * J / (s - discard)))
    refs[[paste(p, n, sep = "_")]] <- list(
      s = s, r = r, spec = spec, pair = pair, chunks = chunks,
      critical = qnorm(1 - cfg$alpha / 2)
    )
  }
  refs
}

# Training cosines are centered within chunks; testing cosines stay raw.
projection_statistics <- function(training, testing, pair, Mvals) {
  Mmax <- max(Mvals)
  w <- scores <- matrix(0, sum(pair), Mmax)
  counts <- numeric(Mmax)
  for (F in training) {
    F <- sweep(F, 2, colMeans(F), "-")
    for (m in seq_len(Mmax)) {
      N <- nrow(F) - m
      w[, m] <- w[, m] + crossprod(
        F[seq_len(N), , drop = FALSE],
        F[seq_len(N) + m, , drop = FALSE]
      )[pair]
      counts[m] <- counts[m] + N
    }
  }
  w <- sweep(w, 2, counts, "/")
  for (m in seq_len(Mmax)) {
    N <- nrow(testing) - m
    scores[, m] <- crossprod(
      testing[seq_len(N), , drop = FALSE],
      testing[seq_len(N) + m, , drop = FALSE]
    )[pair] / sqrt(N)
  }
  vapply(Mvals, function(M) {
    direction <- w[, seq_len(M), drop = FALSE]
    norm <- sqrt(sum(direction^2))
    if (norm == 0) return(scores[1L, 1L])
    sum(direction * scores[, seq_len(M), drop = FALSE]) / norm
  }, numeric(1L))
}

simulate_fitted_replication <- function(p, n, irep, cfg, refs) {
  ref <- refs[[paste(p, n, sep = "_")]]
  blocks <- c(list(seq.int(ref$s + 1L, n)), ref$chunks)
  set_stream_seed(0L, p, n, irep, 3L)
  grids <- lapply(lengths(blocks), function(k)
    matrix(runif(k * p, -1, 1), k, p))
  features <- lapply(grids, function(G)
    sweep(basis_on_grid(G, ref$spec), 2, sqrt(ref$spec$lambda), "/"))
  statistic <- matrix(0, length(cfg$cases), length(cfg$M))
  for (ii in seq_along(cfg$cases)) {
    case <- cfg$cases[ii]
    set_stream_seed(FITTED_STREAM_ID[case], p, n, irep, 2L)
    obj <- generate_case(case, n, p, cfg$burnin, cfg$smooth_exp_tau)
    Z <- fit_residuals(obj, ref$s)
    F <- lapply(seq_along(blocks), function(j) {
      matching <- ot_assignment(Z[blocks[[j]], , drop = FALSE], grids[[j]])
      features[[j]][matching, , drop = FALSE]
    })
    statistic[ii, ] <- projection_statistics(F[-1L], F[[1L]],
                                             ref$pair, cfg$M)
  }
  data.frame(
    case = rep(cfg$cases, each = length(cfg$M)), p = p, n = n,
    s = ref$s, r = ref$r, rep = irep,
    M = rep(cfg$M, times = length(cfg$cases)),
    statistic = as.vector(t(statistic)), critical = ref$critical,
    reject = as.vector(t(abs(statistic) > ref$critical))
  )
}

# -----------------------------------------------------------------------------
# Parallel runs, checkpoints, and summaries
# -----------------------------------------------------------------------------

summarize_results <- function(raw) {
  groups <- interaction(raw$case, raw$p, raw$n, raw$test, raw$M, drop = TRUE)
  ans <- do.call(rbind, lapply(split(raw, groups), function(x) {
    rate <- mean(x$reject)
    cbind(x[1L, c("case", "p", "n", "test", "M")],
          replications = nrow(x), rejection_rate = rate,
          mc_se = sqrt(rate * (1 - rate) / nrow(x)))
  }))
  rownames(ans) <- NULL
  ans[order(ans$case, ans$p, ans$n, ans$test, ans$M), ]
}

run_experiment <- function(cfg, core_file) {
  dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
  refs <- if (cfg$setting == "fitted") prepare_fitted_reference(cfg) else
    prepare_observed_reference(cfg)
  simulate <- switch(cfg$setting, observed = simulate_observed_replication,
                     fitted = simulate_fitted_replication)
  total <- cfg$reps
  # The hash prevents mixing checkpoints from different code or settings.
  stamp <- tempfile()
  saveRDS(list(code = tools::md5sum(core_file),
               cfg = cfg[setdiff(names(cfg), c("ncores", "outdir"))]), stamp)
  checkpoint <- file.path(cfg$outdir, paste0("checkpoints_", tools::md5sum(stamp)))
  unlink(stamp)
  dir.create(checkpoint, showWarnings = FALSE)
  saveRDS(list(cfg = cfg, refs = refs), file.path(cfg$outdir, "settings.rds"))
  cluster <- parallel::makeCluster(cfg$ncores)
  on.exit(parallel::stopCluster(cluster))
  parallel::clusterCall(cluster, function(path, libs) {
    .libPaths(libs)
    source(path)
    NULL
  }, normalizePath(core_file), .libPaths())
  parallel::clusterExport(cluster, c("cfg", "refs", "simulate"), envir = environment())
  design <- expand.grid(p = cfg$p, n = cfg$n)
  blocks <- vector("list", nrow(design))
  timing_file <- file.path(checkpoint, "timings.csv")
  timings <- if (file.exists(timing_file)) read.csv(timing_file) else
    transform(design, seconds = 0)

  for (j in seq_len(nrow(design))) {
    p <- design$p[j]
    n <- design$n[j]
    file <- file.path(checkpoint, sprintf("p%d_n%d.rds", p, n))
    raw <- if (file.exists(file)) readRDS(file) else NULL
    done <- if (is.null(raw)) 0L else max(raw$rep)
    start <- proc.time()[3L]
    while (done < total) {
      ids <- seq.int(done + 1L, min(total, done + cfg$batch_size))
      rows <- parallel::parLapplyLB(cluster, ids, function(i, p, n)
        simulate(p, n, i, cfg, refs), p = p, n = n)
      raw <- rbind(raw, do.call(rbind, rows))
      done <- max(ids)
      saveRDS(raw, file, compress = FALSE)
      message(sprintf("%s p=%d n=%d: %d/%d (%.1f min)",
                      cfg$setting, p, n, done, total,
                      (proc.time()[3L] - start) / 60))
    }
    timings$seconds[j] <- timings$seconds[j] + proc.time()[3L] - start
    write.csv(timings, timing_file, row.names = FALSE)
    blocks[[j]] <- raw
  }

  write.csv(timings, file.path(cfg$outdir, "timings.csv"), row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(cfg$outdir, "sessionInfo.txt"))
  raw <- do.call(rbind, blocks)
  if (cfg$setting == "fitted") raw$test <- "projection"
  summary <- summarize_results(raw)
  saveRDS(raw, file.path(cfg$outdir, "raw_results.rds"))
  write.csv(summary, file.path(cfg$outdir, "summary_results.csv"), row.names = FALSE)
  writexl::write_xlsx(list(Summary = summary, Timings = timings),
                      file.path(cfg$outdir, "simulation_results.xlsx"))
  invisible(summary)
}
