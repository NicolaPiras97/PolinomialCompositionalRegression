if(!require(lpSolveAPI)) install.packages("lpSolveAPI")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(codalm)) install.packages("codalm") 
if(!require(gtools)) install.packages("gtools")  
if(!require(compositions)) install.packages("compositions") 
if(!require(dplyr)) install.packages("dplyr") 
if(!require(tidyr)) install.packages("tidyr") 
if(!require(reshape2)) install.packages("reshape2") 

generate_data_poly <- function(
    N = 200,
    Cx = 5,
    Cy = 8,
    degree = 2,
    response_type = c("ordinal", "nominal"),
    lambda = 1.5,
    alpha_dir = 60,
    B=NULL,
    x_weights = NULL,
    y_scores = NULL,
    beta = NULL,
    normalize_eta = TRUE
) {
  
  
  
  # ============================================================
  # 1. Generazione X composizionale
  # ============================================================
  
  X <- matrix(0, N, Cx)
  
  for (n in 1:N) {
    X[n, ] <- rdirichlet(
      1,
      rep(1, Cx)
    )
  }
  
  
  # ============================================================
  # 2. Costruzione della base polinomiale di Bernstein
  # ============================================================
  
  Z <- bernstein_basis(
    X,
    degree = degree,
    alpha_weights = x_weights
  )
  
  D <- ncol(Z)
  
  if (is.null(B)) {
    B <- matrix(
      rnorm(D * Cy, mean = 0, sd = 0.5),
      nrow = D,
      ncol = Cy
    )
  } else {
    if (!all(dim(B) == c(D, Cy))) {
      stop(
        sprintf(
          "B must be a %d x %d matrix",
          D, Cy
        )
      )
    }
  }
  
  
  # ============================================================
  # 3. Coefficienti veri del polinomio
  # ============================================================
  
  if (is.null(beta)) {
    
    beta <- rnorm(
      D,
      mean = 0,
      sd = 1
    )
    
  } else {
    
    if (length(beta) != D) {
      stop(
        sprintf(
          "beta must have length %d for Cx = %d and degree = %d",
          D, Cx, degree
        )
      )
    }
  }
  
  
  # ============================================================
  # 4. Predittore polinomiale
  # ============================================================
  
  eta_raw <- as.vector(
    Z %*% beta
  )
  
  
  # ============================================================
  # 5. Riscalamento del predittore
  #    nello spazio delle categorie ordinali
  # ============================================================
  
  if (is.null(y_scores)) {
    zy <- 1:Cy
  } else {
    
    if (length(y_scores) != Cy) {
      stop("y_scores must have length Cy")
    }
    
    zy <- y_scores
  }
  
  
  if (normalize_eta) {
    
    eta_min <- min(eta_raw)
    eta_max <- max(eta_raw)
    
    if (eta_max == eta_min) {
      eta <- rep(mean(zy), N)
    } else {
      
      eta <- min(zy) +
        (eta_raw - eta_min) /
        (eta_max - eta_min) *
        (max(zy) - min(zy))
    }
    
  } else {
    
    eta <- eta_raw
  }
  
  
  # ============================================================
  # 6. Generazione Y composizionale
  #    tramite kernel Laplace
  # ============================================================
  
  Y <- matrix(
    0,
    nrow = N,
    ncol = Cy
  )
  
  for (n in 1:N) {
    
    # kernel Laplace
    if (response_type == "ordinal") {
      
      # Ordinal mechanism
      prob <- exp(
        -lambda * abs(zy - eta[n])
      )
      
      # Normalization
      prob <- prob / sum(prob)
      
    } else {
      
      # Non-ordinal mechanism
      theta <- as.vector(Z[n, ] %*% B)
      
      # Numerical stabilization of the softmax
      prob <- exp(theta - max(theta))
      
      # Normalization
      prob <- prob / sum(prob)
    }
    
    # Y composizionale
    Y[n, ] <- rdirichlet(
      1,
      alpha_dir * prob
    )
  }
  
  
  # ============================================================
  # 7. Output
  # ============================================================
  
  return(
    list(
      X = X,
      Y = Y,
      Z = Z,
      beta = beta,
      eta = eta,
      eta_raw = eta_raw,
      y_scores = zy
    )
  )
}

run_simulation <- function(
    iter = 1000,
    N = 100,
    Cx0 = 5,
    Cy0 = 8,
    distance = "wd",
    response_type="ordinal",
    degree = 2,
    lambda = 2,
    alpha_dir = 80,
    train_ratio = 0.7,
    prop_zero_max = 0,
    prop_zero_max2 = 0,
    prop_zero =0,
    prop_inv=0
){
  
  result <- matrix(0, nrow = iter, ncol = 14)
  colnames(result) <- c(
    "RPS_OT_OT","RPS_COD_OT",
    "Brier_OT_OT","Brier_COD_OT",
    "RPS_OT_COD","RPS_COD_COD",
    "Brier_OT_COD","Brier_COD_COD",
    "Win_RPS_OT_DGP","Win_Brier_OT_DGP",
    "Win_RPS_COD_DGP","Win_Brier_COD_DGP",
    "Win_RPS_TRUE_DGP","Win_Brier_TRUE_DGP"
  )
  
  for(it in 1:iter){
    
    # --- GENERAZIONE DGP ---
    data <- generate_data_poly(
      N = N,
      Cx = Cx0,
      Cy = Cy0,
      degree=degree,
      response_type = response_type,
      lambda = lambda,
      alpha_dir = alpha_dir
    )
    
    x0 <- data$Z
    Cx<-ncol(x0)
    y  <- data$Y
    
    # --- TRAIN/VALID SPLIT ---
    tr <- trunc(train_ratio * N)
    va <- N - tr
    
    idx  <- sample(1:N, N)
    intr <- idx[1:tr]
    inva <- idx[(tr+1):N]
    
    # --- FORZO PERCENTUALE DI ZERI NEL VALIDATION SET ---
    if(prop_zero_max > 0){
      
      # numero di righe da modificare
      n_rows <- length(inva)
      n_zero_rows <- floor(prop_zero_max* n_rows)
      
      if(n_zero_rows > 0){
        
        # selezione casuale delle righe nel validation set
        selected_rows <- sample(inva, n_zero_rows, replace = FALSE)
        
        for(r in selected_rows){
          
          # indice della componente più grande
          j_max <- which.max(y[r,])
          
          y[r, j_max] <- 0
          
          if(sum(y[r,]) > 0){
            y[r,] <- y[r,] / sum(y[r,])
          }
        }
      }
    }
    
    if(prop_zero_max2 > 0){
      
      # numero di righe da modificare
      n_rows <- length(inva)
      n_zero_rows <- floor(prop_zero_max2 * n_rows)
      
      if(n_zero_rows > 0){
        
        # selezione casuale delle righe nel validation set
        selected_rows <- sample(inva, n_zero_rows, replace = FALSE)
        
        for(r in selected_rows){
          
          # controllo numero colonne
          if(ncol(y) >= 5){
            # indici delle due componenti più grandi
            j_max <- order(y[r, ], decreasing = TRUE)[1:2]
          } else {
            # indice della componente più grande
            j_max <- which.max(y[r, ])
          }
          
          # azzera le componenti selezionate
          y[r, j_max] <- 0
          
          # normalizzazione se la somma è positiva
          if(sum(y[r,]) > 0){
            y[r,] <- y[r,] / sum(y[r,])
          }
        }
      }
    }
    
    if(prop_zero > 0){
      
      n_rows <- length(inva)
      n_zero_rows <- floor(prop_zero * n_rows)
      
      if(n_zero_rows > 0){
        
        selected_rows <- sample(inva, n_zero_rows, replace = FALSE)
        
        for(r in selected_rows){
          
          n_comp <- ncol(y)
          
          # quante componenti annullare
          k <- if(n_comp >= 5) 2 else 1
          
          # indici candidati (solo quelli > 0 per evitare inutilità)
          non_zero_idx <- which(y[r, ] > 0)
          
          # se ci sono abbastanza componenti non nulle
          if(length(non_zero_idx) >= k){
            
            j_sel <- sample(non_zero_idx, k, replace = FALSE)
            
          } else {
            # fallback: prendi da tutte le colonne
            j_sel <- sample(seq_len(n_comp), k, replace = FALSE)
          }
          
          # annulla le componenti selezionate
          y[r, j_sel] <- 0
          
          # rinormalizza
          s <- sum(y[r,])
          if(s > 0){
            y[r,] <- y[r,] / s
          } else {
            # fallback opzionale (evita riga tutta zero)
            y[r,] <- rep(1/ncol(y), ncol(y))
          }
        }
      }
    }
    
    if(prop_inv > 0){
      
      n_rows <- length(inva)
      n_inv_rows <- floor(prop_inv * n_rows)
      
      if(n_inv_rows > 0){
        
        selected_rows <- sample(inva, n_inv_rows, replace = FALSE)
        
        for(r in selected_rows){
          
          # inversione delle componenti
          y[r, ] <- rev(y[r, ])
        }
      }
    }
    
    # COD matrix from DGP
    Ac <- codalm(y, x0)
    
    weights <- rep(1, Cy0-1)
    
    # OT matrix from DGP
    ydata <- lapply(1:N, function(i) y[i,])
    xdata <- lapply(1:N, function(i) x0[i,])
    if(distance=="wd"){
    sol <- solve_simplex_lp(xdata, ydata, weights)
    } else if(distance=="cvm"){
    sol <- fit_simplex_ordinal(x0,y,degree=1,distance="cvm")
    } else if(distance=="scs"){
    sol <- fit_simplex_nominal(x0,y,degree=1)
    } else {
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    mat   <- sol$A
    
    # --- GENERO NUOVO CAMPIONE ---
    x  <- matrix(0,N,Cx)
    y1 <- matrix(0,N,Cy0)
    y2 <- matrix(0,N,Cy0)
    
    for(n in 1:N){
      x[n,]  <- rdirichlet(1, rep(1,Cx))
      y1[n,] <- safe_dirichlet(mat %*% x[n,])
      y2[n,] <- safe_dirichlet(t(Ac) %*% x[n,])
    }
    
    
    Plist  <- lapply(1:N, function(i) x[i,])
    P1list <- lapply(1:N, function(i) y1[i,])
    P2list <- lapply(1:N, function(i) y2[i,])
    
    Plist_tr  <- lapply(intr, function(i) x[i,])
    P1list_tr <- lapply(intr, function(i) y1[i,])
    P2list_tr <- lapply(intr, function(i) y2[i,])
    
    # --- STIMA SU TRAIN ---
    if(distance=="wd"){
      A1 <- solve_simplex_lp(Plist_tr, P1list_tr, weights)$A
    } else if(distance=="cvm"){
      A1 <- fit_simplex_ordinal(x,y1,degree=1,distance="cvm")$A
    } else if(distance=="scs"){
      A1 <- fit_simplex_nominal(x,y1,degree=1)$A
    } else {
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    if(distance=="wd"){
      A2 <- solve_simplex_lp(Plist_tr, P2list_tr, weights)$A
    } else if(distance=="cvm"){
      A2 <- fit_simplex_ordinal(x,y2,degree=1,distance="cvm")$A
    } else if(distance=="scs"){
      A2 <- fit_simplex_nominal(x,y2,degree=1)$A
    } else {
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    # COD stimato
    eps <- 1e-8
    Ymat1 <- do.call(rbind, P1list_tr)
    Ymat2 <- do.call(rbind, P2list_tr)
    Xmat <- do.call(rbind, Plist_tr)
    Ymat1[Ymat1 < eps] <- eps
    Ymat2[Ymat2 < eps] <- eps
    Xmat[Xmat < eps] <- eps
    Ymat1 <- Ymat1 / rowSums(Ymat1)
    Ymat2 <- Ymat2 / rowSums(Ymat2)
    Xmat <- Xmat / rowSums(Xmat)
    A1_cod <- codalm(Ymat1, Xmat)
    A2_cod <- codalm(Ymat2, Xmat)
    
    # --- ERRORI ---
    calc_errors <- function(A_ot, A_cod, y_true, x_mat){
      err_wrps_ot <- err_wrps_cod <- err_b_ot <- err_b_cod <- err_w1_ot <- err_w1_cod <- 0
      for(i in inva){
        pred_ot  <- as.vector(A_ot  %*% x_mat[i,])
        pred_cod <- as.vector(t(A_cod) %*% x_mat[i,])
        err_wrps_ot  <- err_wrps_ot  + calc_wrps(pred_ot,  y_true[i,])
        err_wrps_cod <- err_wrps_cod + calc_wrps(pred_cod, y_true[i,])
        err_b_ot     <- err_b_ot  + calc_brier(pred_ot,  y_true[i,])
        err_b_cod    <- err_b_cod + calc_brier(pred_cod, y_true[i,])
      }
      return(c(
        err_wrps_ot/va, err_wrps_cod/va,
        err_b_ot/va,    err_b_cod/va
      ))
    }
    
    errs1 <- calc_errors(A1, A2_cod, y1, x)
    errs2 <- calc_errors(A2, A1_cod, y2, x)
    result[it,1:4]  <- errs1
    result[it,5:8] <- errs2
    
    # --- WINS ---
    result[it,9:10] <- as.numeric(errs1[c(1,3)] < errs1[c(2,4)])
    result[it,11:12] <- as.numeric(errs2[c(1,3)] < errs2[c(2,4)])
    result[it,13:14] <- as.numeric(errs1[c(1,3)] < errs2[c(2,4)])
  }
  
  list(
    mean_res  = colMeans(result),
    sd_res    = apply(result,2,sd),
    win_rates = colMeans(result[,9:14]),
    raw_results = result
  )
}

s1 <- run_simulation(
  iter=10,
  N = 50,
  Cx0 = 5,
  Cy0 = 5,
  response_type="nominal",
  distance="scs"
)
s1$win_rates         

###########################################################################
run_simulation_caseI <- function(
    iter = 1000,
    N = 100,
    Cx0 = 5,
    Cy0 = 8,
    response_type = "ordinal",
    degree = 2,
    lambda = 2,
    alpha_dir = 80,
    tol_equal = 1e-8
){
  
  result <- matrix(NA, nrow = iter, ncol = 6)
  
  colnames(result) <- c(
    "W1",
    "SCS",
    "CVM",
    "R2_W1",
    "R2_SCS",
    "R2_CVM"
    #"OPI_W1",
    #"OPI_CVM"
  )
  
  for(it in 1:iter){
    
    # =====================================================
    # GENERAZIONE DGP
    # =====================================================
    
    data <- generate_data_poly(
      N = N,
      Cx = Cx0,
      Cy = Cy0,
      degree = degree,
      response_type = response_type,
      lambda = lambda,
      alpha_dir = alpha_dir
    )
    
    x0 <- data$Z
    Cx <- ncol(x0)
    y  <- data$Y
    
    X <- x0
    Y <- y
    
    
    # =====================================================
    # WEIGHTS
    # =====================================================
    
    weights_unit <- rep(1, Cy0 - 1)
    
    wsum_u <- sum(weights_unit)
    wsum_cvm <- sum(weights_unit^2)
    
    # =====================================================
    # LISTE
    # =====================================================
    
    X_list <- lapply(1:N, function(i) X[i,])
    Y_list <- lapply(1:N, function(i) Y[i,])
    
    
    # =====================================================
    # W1
    # =====================================================
    
    A_unit <- solve_simplex_lp(
      X_list,
      Y_list,
      weights_unit
    )$A
    
    
    SSE_u <- 0
    
    pred_u <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_u[i,] <- as.vector(
        A_unit %*% X[i,]
      )
    }
    
    for(i in 1:N){
      SSE_u <- SSE_u +
        wd(
          weights_unit,
          Y[i,],
          pred_u[i,]
        )
    }
    
    SSE_u <- SSE_u
    
    result[it,"W1"] <- SSE_u/(wsum_u*N)
    
    
    result[it,"R2_W1"] <- compute_R2(
      Y,
      X,
      A_unit,
      distance = "wd",
      weights = weights_unit
    )$R2
    
    
    # =====================================================
    # SCS
    # =====================================================
    
    A_scs <- fit_simplex_nominal(
      X,
      Y,
      degree = 1
    )$A
    
    
    pred_scs <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_scs[i,] <- as.vector(
        A_scs %*% X[i,]
      )
    }
    
    SSE_scs <- 0
    
    for(i in 1:N){
      SSE_scs <- SSE_scs +
        sym_chi2(
          Y[i,],
          pred_scs[i,]
        )
    }
    
    result[it,"SCS"] <- SSE_scs/(4*N)
    
    
    result[it,"R2_SCS"] <- compute_R2(
      Y,
      X,
      A_scs,
      distance = "scs"
    )$R2
    
    
    # =====================================================
    # CVM
    # =====================================================
    
    A_cvm <- fit_simplex_ordinal(
      X,
      Y,
      degree = 1,
      distance = "cvm"
    )$A
    
    
    pred_cvm <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_cvm[i,] <- as.vector(
        A_cvm %*% X[i,]
      )
    }
    
    SSE_cvm <- 0
    
    for(i in 1:N){
      SSE_cvm <- SSE_cvm +
        cvm_distance(
          weights_unit,
          Y[i,],
          pred_cvm[i,]
        )
    }
    
    result[it,"CVM"] <- SSE_cvm/(wsum_cvm*N)
    
    
    result[it,"R2_CVM"] <- compute_R2(
      Y,
      X,
      A_cvm,
      distance = "cvm",
      weights = weights_unit
    )$R2
    
    
    # =====================================================
    # OPI W1
    # =====================================================
    
    #ind_u <- 0
    
    #for(i in 1:(N-1)){
      
    #  for(j in (i+1):N){
        
    #    if(
    #      opiwd(
    #        weights_unit,
    #        Y[i,],
    #        Y[j,]
    #      ) ==
    #      opiwd(
    #        weights_unit,
    #        pred_u[i,],
    #        pred_u[j,]
    #      )
    #    ){
    #      ind_u <- ind_u + 1
    #    }
    #  }
    #}
    
    #result[it,"OPI_W1"] <-
    #  ind_u / choose(N,2)
    
    
    # =====================================================
    # OPI CVM
    # =====================================================
    
    #ind_cvm <- 0
    
    #for(i in 1:(N-1)){
      
    #  for(j in (i+1):N){
        
    #    d_obs <- cvm_distance(
    #      weights_unit,
    #      Y[i,],
    #      Y[j,]
    #    )
        
    #    d_pred <- cvm_distance(
    #      weights_unit,
    #      pred_cvm[i,],
    #      pred_cvm[j,]
    #    )
        
    #    if(
    #      abs(d_obs - d_pred) <= tol_equal
    #    ){
    #      ind_cvm <- ind_cvm + 1
    #    }
    #  }
    #}
    
    #result[it,"OPI_CVM"] <-
    #  ind_cvm / choose(N,2)
  }
  
  
  list(
    mean_res = colMeans(
      result,
      na.rm = TRUE
    ),
    
    sd_res = apply(
      result,
      2,
      sd,
      na.rm = TRUE
    ),
    
    raw_results = result
  )
}

s1 <- run_simulation_caseI(
  iter=10,
  N = 100,
  Cx = 3,
  Cy = 3,
  response_type = "nominal"
)

