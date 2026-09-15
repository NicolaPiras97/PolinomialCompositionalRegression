safe_dirichlet <- function(alpha, scale = 10, eps = 1e-6){
  alpha <- scale * alpha
  alpha[alpha < eps] <- eps
  alpha[!is.finite(alpha)] <- eps
  return(rdirichlet(1, alpha))
}

calc_wrps <- function(p_hat, y_true){
  Fhat <- cumsum(p_hat)
  Fobs <- cumsum(y_true)
  sum((Fhat - Fobs)^2)
}

# Brier multicategoria
calc_brier <- function(p_hat, y_true){
  sum((p_hat - y_true)^2)
}

# Wasserstein W1 per ordinale
calc_wass <- function(p_hat, y_true){
  Fhat <- cumsum(p_hat)
  Fobs <- cumsum(y_true)
  sum(abs(Fhat - Fobs))
}   

wd <- function(weights, P, Q){
  n <- length(P)
  sum(weights * abs(cumsum(P)[1:(n-1)] - cumsum(Q)[1:(n-1)]))
}


# ---------------------------------------------------------
# Symmetric Chi-squared distance
# ---------------------------------------------------------
sym_chi2 <- function(P, Q) {
  
  if (length(P) != length(Q)) {
    stop("P and Q must have the same length.")
  }
  
  denominator <- (P + Q) / 2
  
  # Avoid division by zero when P_j = Q_j = 0
  idx <- denominator > 0
  
  sum((P[idx] - Q[idx])^2 / denominator[idx])
}


# ---------------------------------------------------------
# Cramér-von Mises distance
# ---------------------------------------------------------
cvm_distance <- function(weights, P, Q) {
  
  if (length(P) != length(Q)) {
    stop("P and Q must have the same length.")
  }
  
  if (length(weights) != length(P) - 1) {
    stop("For the Cramér-von Mises distance, weights must have length m.")
  }
  
  # Cumulative distributions, excluding the last category
  FP <- cumsum(P)[1:(length(P) - 1)]
  FQ <- cumsum(Q)[1:(length(Q) - 1)]
  
  sum(weights^2 * (FP - FQ)^2)
}


# ---------------------------------------------------------
# General distance function
# ---------------------------------------------------------
compute_distance <- function(P, Q, distance = c("wd",
                                                "scs",
                                                "cvm"),
                             weights = NULL) {
  
  distance <- match.arg(distance)
  
  if (distance == "wd") {
    
    if (is.null(weights)) {
      stop("weights must be provided for the Wasserstein distance.")
    }
    
    return(wd(weights, P, Q))
    
  } else if (distance == "scs") {
    
    return(sym_chi2(P, Q))
    
  } else if (distance == "cvm") {
    
    if (is.null(weights)) {
      stop("weights must be provided for the Cramér-von Mises distance.")
    }
    
    return(cvm_distance(weights, P, Q))
  }
}


# ---------------------------------------------------------
# General R2
# ---------------------------------------------------------
compute_R2 <- function(Y, X, A,
                       distance = c("wd",
                                    "scs",
                                    "cvm"),
                       weights = NULL) {
  
  distance <- match.arg(distance)
  
  N <- nrow(Y)
  Cy <- ncol(Y)
  
  # Wasserstein-Fréchet mean
  ymean <- compute_wfrechet_mean(Y)
  
  # Predicted responses
  pred <- matrix(0, nrow = N, ncol = Cy)
  
  for (i in 1:N) {
    pred[i, ] <- as.vector(A %*% X[i, ])
  }
  
  # Residual and total distances
  Dres <- 0
  Dtot <- 0
  
  for (i in 1:N) {
    
    Dres <- Dres + compute_distance(
      Y[i, ],
      pred[i, ],
      distance = distance,
      weights = weights
    )
    
    Dtot <- Dtot + compute_distance(
      Y[i, ],
      ymean,
      distance = distance,
      weights = weights
    )
  }
  
  if (Dtot == 0) {
    return(NA)
  }
  
  R2 <- 1 - Dres / Dtot
  
  return(list(
    R2 = R2,
    Dres = Dres,
    Dtot = Dtot,
    y_mean = ymean,
    pred = pred
  ))
}
