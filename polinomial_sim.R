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
    x_weights = NULL,
    y_scores = NULL,
    beta = NULL,
    B = NULL
) {
  
  response_type <- match.arg(response_type)
  
  # ---------------------------------------------------------
  # 1. Generazione di X composizionale
  # ---------------------------------------------------------
  X <- matrix(0, N, Cx)
  
  for (n in 1:N) {
    X[n, ] <- rdirichlet(1, rep(1, Cx))
  }
  
  # ---------------------------------------------------------
  # 2. Base polinomiale di Bernstein
  # ---------------------------------------------------------
  Z <- bernstein_basis(
    X,
    degree = degree,
    alpha_weights = x_weights
  )
  
  D <- ncol(Z)
  
  # =========================================================
  # CASO ORDINALE
  # =========================================================
  if (response_type == "ordinal") {
    
    # -------------------------------------------------------
    # 3. Coefficienti del predittore latente
    #
    #    eta_i = Z_i %*% beta
    #
    #    Di default:
    #    beta = 1, 2, ..., D
    # -------------------------------------------------------
    if (is.null(beta)) {
      beta <- 1:D
    } else {
      if (length(beta) != D) {
        stop(sprintf(
          "beta must have length %d for Cx = %d and degree = %d",
          D, Cx, degree
        ))
      }
    }
    
    # Predittore latente polinomiale
    eta <- as.vector(Z %*% beta)
    
    # -------------------------------------------------------
    # 4. Coordinate latenti delle categorie ordinali di Y
    # -------------------------------------------------------
    if (is.null(y_scores)) {
      zy <- 1:Cy
    } else {
      if (length(y_scores) != Cy) {
        stop("y_scores must have length Cy")
      }
      zy <- y_scores
    }
    
    # -------------------------------------------------------
    # 5. Generazione di Y | X
    #
    #    p_ic propto exp(-lambda * |zy_c - eta_i|)
    # -------------------------------------------------------
    Y <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      prob <- exp(-lambda * abs(zy - eta[n]))
      prob <- prob / sum(prob)
      
      Y[n, ] <- rdirichlet(
        1,
        alpha_dir * prob
      )
    }
    
    return(list(
      X = X,
      Y = Y,
      Z = Z,
      beta = beta,
      eta = eta,
      y_scores = zy
    ))
  }
  
  
  # =========================================================
  # CASO NOMINALE
  # =========================================================
  if (response_type == "nominal") {
    
    # -------------------------------------------------------
    # 3. Matrice dei coefficienti polinomiali
    #
    #    B ha dimensione D x Cy.
    #
    #    Ogni colonna corrisponde a una categoria di Y.
    # -------------------------------------------------------
    if (is.null(B)) {
      
      B <- matrix(
        rnorm(D * Cy, mean = 0, sd = 0.5),
        nrow = D,
        ncol = Cy
      )
      
    } else {
      
      if (!all(dim(B) == c(D, Cy))) {
        stop(sprintf(
          "B must be a %d x %d matrix",
          D, Cy
        ))
      }
    }
    
    # -------------------------------------------------------
    # 4. Predittori polinomiali per le Cy categorie
    #
    #    theta_ic = Z_i %*% B[, c]
    # -------------------------------------------------------
    theta <- Z %*% B
    
    # -------------------------------------------------------
    # 5. Softmax + generazione composizionale di Y | X
    # -------------------------------------------------------
    Y <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      # Softmax con stabilizzazione numerica
      prob <- exp(theta[n, ] - max(theta[n, ]))
      prob <- prob / sum(prob)
      
      Y[n, ] <- rdirichlet(
        1,
        alpha_dir * prob
      )
    }
    
    return(list(
      X = X,
      Y = Y,
      Z = Z,
      B = B,
      theta = theta
    ))
  }
}

generate_linear_poly_data <- function(
    N = 200,
    Cx = 5,
    Cy = 8,
    degree = 2,
    response_type = c("ordinal", "nominal"),
    lambda = 1.5,
    alpha_dir = 60,
    x_weights = NULL,
    y_scores = NULL,
    beta = NULL,
    B = NULL,
    sigma_B = 0.5
) {
  
  response_type <- match.arg(response_type)
  
  # ============================================================
  # 1. GENERATE X
  # ============================================================
  
  X <- matrix(0, N, Cx)
  
  for (n in 1:N) {
    X[n, ] <- rdirichlet(1, rep(1, Cx))
  }
  
  
  # ============================================================
  # 2. BERNSTEIN BASIS
  # ============================================================
  
  Z <- bernstein_basis(
    X,
    degree = degree,
    alpha_weights = x_weights
  )
  
  D <- ncol(Z)
  
  
  # ============================================================
  # 3. LINEAR + POLYNOMIAL
  # ============================================================
  
  if (response_type == "ordinal") {
    
    # ------------------------------------------------------------
    # ORDINAL CASE
    # ------------------------------------------------------------
    
    # Linear latent score
    zx <- 1:Cx
    s <- as.vector(X %*% zx)
    
    
    # Scores of Y categories
    if (is.null(y_scores)) {
      zy <- 1:Cy
    } else {
      if (length(y_scores) != Cy) {
        stop("y_scores must have length Cy")
      }
      zy <- y_scores
    }
    
    
    # Polynomial latent score
    if (is.null(beta)) {
      beta <- 1:D
    } else {
      if (length(beta) != D) {
        stop(sprintf(
          "beta must have length %d for Cx = %d and degree = %d",
          D, Cx, degree
        ))
      }
    }
    
    eta <- as.vector(Z %*% beta)
    
    
    # ------------------------------------------------------------
    # Generate Y: linear
    # ------------------------------------------------------------
    
    Y_linear <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      prob <- exp(
        -lambda * abs(zy - s[n])
      )
      
      prob <- prob / sum(prob)
      
      Y_linear[n, ] <- rdirichlet(
        1,
        alpha_dir * prob
      )
    }
    
    
    # ------------------------------------------------------------
    # Generate Y: polynomial
    # ------------------------------------------------------------
    
    Y_poly <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      prob <- exp(
        -lambda * abs(zy - eta[n])
      )
      
      prob <- prob / sum(prob)
      
      Y_poly[n, ] <- rdirichlet(
        1,
        alpha_dir * prob
      )
    }
    
    
    # ------------------------------------------------------------
    # Return
    # ------------------------------------------------------------
    
    return(list(
      X = X,
      
      # Bernstein basis
      Z = Z,
      
      # Linear model
      s = s,
      Y_linear = Y_linear,
      
      # Polynomial model
      beta = beta,
      eta = eta,
      Y_poly = Y_poly,
      
      # Y latent category scores
      y_scores = zy
    ))
  }
  
  
  # ============================================================
  # NOMINAL CASE
  # ============================================================
  
  if (response_type == "nominal") {
    
    # ------------------------------------------------------------
    # LINEAR CASE
    # ------------------------------------------------------------
    
    # Coefficient matrix for the linear model
    #
    # Each column corresponds to one Y category.
    #
    # B_linear[j,c] = effect of X_j on category c
    #
    
    if (is.null(B)) {
      
      B_linear <- matrix(
        rnorm(Cx * Cy, mean = 0, sd = sigma_B),
        nrow = Cx,
        ncol = Cy
      )
      
    } else {
      
      if (!is.matrix(B) ||
          nrow(B) != Cx ||
          ncol(B) != Cy) {
        
        stop(sprintf(
          "B must be a %d x %d matrix for Cx = %d and Cy = %d",
          Cx, Cy, Cx, Cy
        ))
      }
      
      B_linear <- B
    }
    
    
    # Linear logits
    theta_linear <- X %*% B_linear
    
    
    # Softmax
    P_linear <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      z <- theta_linear[n, ]
      
      # numerical stability
      z <- z - max(z)
      
      prob <- exp(z)
      prob <- prob / sum(prob)
      
      P_linear[n, ] <- prob
    }
    
    
    # Generate compositional Y
    Y_linear <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      Y_linear[n, ] <- rdirichlet(
        1,
        alpha_dir * P_linear[n, ]
      )
    }
    
    
    # ------------------------------------------------------------
    # POLYNOMIAL CASE
    # ------------------------------------------------------------
    
    # B_poly has one coefficient for every
    # Bernstein basis function and every Y category.
    #
    # Dimension: D x Cy
    
    if (is.null(beta)) {
      B_poly <- matrix(
        rnorm(D * Cy, mean = 0, sd = sigma_B),
        nrow = D,
        ncol = Cy
      )
    } else {
      
      # For nominal we use beta as B_poly
      if (!is.matrix(beta) ||
          nrow(beta) != D ||
          ncol(beta) != Cy) {
        
        stop(sprintf(
          "For response_type = 'nominal', beta must be a %d x %d matrix",
          D, Cy
        ))
      }
      
      B_poly <- beta
    }
    
    
    # Polynomial logits
    theta_poly <- Z %*% B_poly
    
    
    # Softmax
    P_poly <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      z <- theta_poly[n, ]
      
      # numerical stability
      z <- z - max(z)
      
      prob <- exp(z)
      prob <- prob / sum(prob)
      
      P_poly[n, ] <- prob
    }
    
    
    # Generate compositional Y
    Y_poly <- matrix(0, N, Cy)
    
    for (n in 1:N) {
      
      Y_poly[n, ] <- rdirichlet(
        1,
        alpha_dir * P_poly[n, ]
      )
    }
    
    
    # ------------------------------------------------------------
    # Return
    # ------------------------------------------------------------
    
    return(list(
      X = X,
      
      # Bernstein basis
      Z = Z,
      
      # Linear model
      B_linear = B_linear,
      theta_linear = theta_linear,
      P_linear = P_linear,
      Y_linear = Y_linear,
      
      # Polynomial model
      B_poly = B_poly,
      theta_poly = theta_poly,
      P_poly = P_poly,
      Y_poly = Y_poly
    ))
  }
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

run_simulation2 <- function(
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
    "RPS_POLY_POLY","RPS_LIN_POLY",
    "Brier_POLY_POLY","Brier_LIN_POLY",
    "RPS_POLY_LIN","RPS_LIN_LIN",
    "Brier_POLY_LIN","Brier_LIN_LIN",
    "Win_RPS_POLY_DGP","Win_Brier_POLY_DGP",
    "Win_RPS_LIN_DGP","Win_Brier_LIN_DGP",
    "Win_RPS_TRUE_DGP","Win_Brier_TRUE_DGP"
  )
  
  for(it in 1:iter){
    
    # --- GENERAZIONE DGP ---
    data <- generate_linear_poly_data(
      N = N,
      Cx = Cx0,
      Cy = Cy0,
      degree=degree,
      response_type = response_type,
      lambda = lambda,
      alpha_dir = alpha_dir
    )
    
    x0_l <- data$X
    Cx_l<-ncol(x0_l)
    y_l  <- data$Y_linear
    x0_p <- data$Z
    Cx_p<-ncol(x0_p)
    y_p  <- data$Y_poly
    
    # --- TRAIN/VALID SPLIT ---
    tr <- trunc(train_ratio * N)
    va <- N - tr
    
    idx  <- sample(1:N, N)
    intr <- idx[1:tr]
    inva <- idx[(tr+1):N]
    
    # ============================================================
    # FORZO PERCENTUALE DI ZERI / INVERSIONI NEL VALIDATION SET
    # ============================================================
    
    # ------------------------------------------------------------
    # Funzione ausiliaria per modificare le righe
    # ------------------------------------------------------------
    
    apply_zero_max <- function(y, rows) {
      
      for(r in rows){
        
        # indice della componente più grande
        j_max <- which.max(y[r, ])
        
        # azzera la componente
        y[r, j_max] <- 0
        
        # rinormalizza
        s <- sum(y[r, ])
        
        if(s > 0){
          y[r, ] <- y[r, ] / s
        } else {
          y[r, ] <- rep(1 / ncol(y), ncol(y))
        }
      }
      
      return(y)
    }
    
    
    apply_zero_max2 <- function(y, rows) {
      
      for(r in rows){
        
        # numero di componenti da azzerare
        if(ncol(y) >= 5){
          j_max <- order(
            y[r, ],
            decreasing = TRUE
          )[1:2]
        } else {
          j_max <- which.max(y[r, ])
        }
        
        # azzera le componenti selezionate
        y[r, j_max] <- 0
        
        # rinormalizza
        s <- sum(y[r, ])
        
        if(s > 0){
          y[r, ] <- y[r, ] / s
        } else {
          y[r, ] <- rep(1 / ncol(y), ncol(y))
        }
      }
      
      return(y)
    }
    
    
    apply_zero_random <- function(y, rows) {
      
      for(r in rows){
        
        n_comp <- ncol(y)
        
        # quante componenti annullare
        k <- if(n_comp >= 5) 2 else 1
        
        # componenti non nulle
        non_zero_idx <- which(y[r, ] > 0)
        
        if(length(non_zero_idx) >= k){
          
          j_sel <- sample(
            non_zero_idx,
            k,
            replace = FALSE
          )
          
        } else {
          
          j_sel <- sample(
            seq_len(n_comp),
            k,
            replace = FALSE
          )
        }
        
        # azzera
        y[r, j_sel] <- 0
        
        # rinormalizza
        s <- sum(y[r, ])
        
        if(s > 0){
          y[r, ] <- y[r, ] / s
        } else {
          y[r, ] <- rep(1 / ncol(y), ncol(y))
        }
      }
      
      return(y)
    }
    
    
    apply_inversion <- function(y, rows) {
      
      for(r in rows){
        y[r, ] <- rev(y[r, ])
      }
      
      return(y)
    }
    
    
    # ============================================================
    # STESSI SOGGETTI PER DGP LINEARE E POLINOMIALE
    # ============================================================
    
    n_rows <- length(inva)
    
    
    # ------------------------------------------------------------
    # 1. Azzera la componente massima
    # ------------------------------------------------------------
    
    if(prop_zero_max > 0){
      
      n_zero_rows <- floor(
        prop_zero_max * n_rows
      )
      
      if(n_zero_rows > 0){
        
        selected_rows <- sample(
          inva,
          n_zero_rows,
          replace = FALSE
        )
        
        # Stesse righe nei due DGP
        y_l <- apply_zero_max(
          y_l,
          selected_rows
        )
        
        y_p <- apply_zero_max(
          y_p,
          selected_rows
        )
      }
    }
    
    
    # ------------------------------------------------------------
    # 2. Azzera le due componenti massime
    # ------------------------------------------------------------
    
    if(prop_zero_max2 > 0){
      
      n_zero_rows <- floor(
        prop_zero_max2 * n_rows
      )
      
      if(n_zero_rows > 0){
        
        selected_rows <- sample(
          inva,
          n_zero_rows,
          replace = FALSE
        )
        
        # Stesse righe nei due DGP
        y_l <- apply_zero_max2(
          y_l,
          selected_rows
        )
        
        y_p <- apply_zero_max2(
          y_p,
          selected_rows
        )
      }
    }
    
    
    # ------------------------------------------------------------
    # 3. Azzera componenti casuali
    # ------------------------------------------------------------
    
    if(prop_zero > 0){
      
      n_zero_rows <- floor(
        prop_zero * n_rows
      )
      
      if(n_zero_rows > 0){
        
        selected_rows <- sample(
          inva,
          n_zero_rows,
          replace = FALSE
        )
        
        # Stesse righe nei due DGP
        y_l <- apply_zero_random(
          y_l,
          selected_rows
        )
        
        y_p <- apply_zero_random(
          y_p,
          selected_rows
        )
      }
    }
    
    
    # ------------------------------------------------------------
    # 4. Inversione delle componenti
    # ------------------------------------------------------------
    
    if(prop_inv > 0){
      
      n_inv_rows <- floor(
        prop_inv * n_rows
      )
      
      if(n_inv_rows > 0){
        
        selected_rows <- sample(
          inva,
          n_inv_rows,
          replace = FALSE
        )
        
        # Stesse righe nei due DGP
        y_l <- apply_inversion(
          y_l,
          selected_rows
        )
        
        y_p <- apply_inversion(
          y_p,
          selected_rows
        )
      }
    }
    
    weights <- rep(1, Cy0-1)
    # LIN matrix from DGP
    ydatac <- lapply(1:N, function(i) y_l[i,])
    xdatac <- lapply(1:N, function(i) x0_l[i,])
    if(distance=="wd"){
      sol <- solve_simplex_lp(xdatac, ydatac, weights)
    } else if(distance=="cvm"){
      sol <- fit_simplex_ordinal(x0_l,y_l,degree=1,distance="cvm")
    } else if(distance=="scs"){
      sol <- fit_simplex_nominal(x0_l,y_l,degree=1)
    } else {
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    Ac   <- sol$A
    
    
    # POLY matrix from DGP
    ydata <- lapply(1:N, function(i) y_p[i,])
    xdata <- lapply(1:N, function(i) x0_p[i,])
    if(distance=="wd"){
      sol <- solve_simplex_lp(xdata, ydata, weights)
    } else if(distance=="cvm"){
      sol <- fit_simplex_ordinal(x0_p,y_p,degree=1,distance="cvm")
    } else if(distance=="scs"){
      sol <- fit_simplex_nominal(x0_p,y_p,degree=1)
    } else {
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    mat   <- sol$A
    
    # --- GENERO NUOVO CAMPIONE ---
    x_l <- matrix(0, N, Cx_l)
    
    for(n in 1:N){
      x_l[n, ] <- rdirichlet(1, rep(1, Cx_l))
    }
    
    # Stessa X, trasformata nella base di Bernstein
    x_p <- bernstein_basis(
      x_l,
      degree = degree
    )
    
    y1 <- matrix(0, N, Cy0)
    y2 <- matrix(0, N, Cy0)
    
    for(n in 1:N){
      
      # POLY
      y1[n, ] <- safe_dirichlet(
        mat %*% x_p[n, ]
      )
      
      # LINEAR / COD
      y2[n, ] <- safe_dirichlet(
        Ac %*% x_l[n, ]
      )
    }
    
    
    # ============================================================
    # LISTE PER LA STIMA
    # ============================================================
    
    # Spazio originale X
    Plist_l <- lapply(1:N, function(i) x_l[i, ])
    
    # Spazio Bernstein Z(X)
    Plist_p <- lapply(1:N, function(i) x_p[i, ])
    
    # Risposte
    P1list <- lapply(1:N, function(i) y1[i, ])  # DGP POLY
    P2list <- lapply(1:N, function(i) y2[i, ])  # DGP LINEARE
    
    
    # Training set
    Plist_l_tr <- lapply(intr, function(i) x_l[i, ])
    Plist_p_tr <- lapply(intr, function(i) x_p[i, ])
    
    P1list_tr <- lapply(intr, function(i) y1[i, ])
    P2list_tr <- lapply(intr, function(i) y2[i, ])
    
    
    # ============================================================
    # STIMA SUL DGP POLINOMIALE
    # ============================================================
    
    # Metodo POLYNOMIAL
    if(distance == "wd"){
      
      A1_poly <- solve_simplex_lp(
        Plist_p_tr,
        P1list_tr,
        weights
      )$A
      
    } else if(distance == "cvm"){
      
      A1_poly <- fit_simplex_ordinal(
        x_p[intr, ],
        y1[intr, ],
        degree = 1,
        distance = "cvm"
      )$A
      
    } else if(distance == "scs"){
      
      A1_poly <- fit_simplex_nominal(
        x_p[intr, ],
        y1[intr, ],
        degree = 1
      )$A
      
    } else {
      
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    
    # Metodo LINEAR
    if(distance == "wd"){
      
      A1_lin <- solve_simplex_lp(
        Plist_l_tr,
        P1list_tr,
        weights
      )$A
      
    } else if(distance == "cvm"){
      
      A1_lin <- fit_simplex_ordinal(
        x_l[intr, ],
        y1[intr, ],
        degree = 1,
        distance = "cvm"
      )$A
      
    } else if(distance == "scs"){
      
      A1_lin <- fit_simplex_nominal(
        x_l[intr, ],
        y1[intr, ],
        degree = 1
      )$A
      
    } else {
      
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    
    # ============================================================
    # STIMA SUL DGP LINEARE
    # ============================================================
    
    # Metodo POLYNOMIAL
    if(distance == "wd"){
      
      A2_poly <- solve_simplex_lp(
        Plist_p_tr,
        P2list_tr,
        weights
      )$A
      
    } else if(distance == "cvm"){
      
      A2_poly <- fit_simplex_ordinal(
        x_p[intr, ],
        y2[intr, ],
        degree = 1,
        distance = "cvm"
      )$A
      
    } else if(distance == "scs"){
      
      A2_poly <- fit_simplex_nominal(
        x_p[intr, ],
        y2[intr, ],
        degree = 1
      )$A
      
    } else {
      
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    
    # Metodo LINEAR
    if(distance == "wd"){
      
      A2_lin <- solve_simplex_lp(
        Plist_l_tr,
        P2list_tr,
        weights
      )$A
      
    } else if(distance == "cvm"){
      
      A2_lin <- fit_simplex_ordinal(
        x_l[intr, ],
        y2[intr, ],
        degree = 1,
        distance = "cvm"
      )$A
      
    } else if(distance == "scs"){
      
      A2_lin <- fit_simplex_nominal(
        x_l[intr, ],
        y2[intr, ],
        degree = 1
      )$A
      
    } else {
      
      stop("distance must be 'wd', 'cvm', or 'scs'")
    }
    
    
    # ============================================================
    # ERRORI
    # ============================================================
    
    calc_errors <- function(
    A_poly,
    A_lin,
    y_true,
    x_poly,
    x_lin
    ){
      
      err_wrps_poly <- 0
      err_wrps_lin  <- 0
      err_b_poly    <- 0
      err_b_lin     <- 0
      
      for(i in inva){
        
        # --------------------------------------------------------
        # Previsione POLY
        # --------------------------------------------------------
        
        pred_poly <- as.vector(
          A_poly %*% x_poly[i, ]
        )
        
        
        # --------------------------------------------------------
        # Previsione LINEARE
        # --------------------------------------------------------
        
        pred_lin <- as.vector(
          A_lin %*% x_lin[i, ]
        )
        
        
        # --------------------------------------------------------
        # RPS
        # --------------------------------------------------------
        
        err_wrps_poly <- err_wrps_poly +
          calc_wrps(
            pred_poly,
            y_true[i, ]
          )
        
        err_wrps_lin <- err_wrps_lin +
          calc_wrps(
            pred_lin,
            y_true[i, ]
          )
        
        
        # --------------------------------------------------------
        # Brier
        # --------------------------------------------------------
        
        err_b_poly <- err_b_poly +
          calc_brier(
            pred_poly,
            y_true[i, ]
          )
        
        err_b_lin <- err_b_lin +
          calc_brier(
            pred_lin,
            y_true[i, ]
          )
      }
      
      
      return(c(
        err_wrps_poly / va,
        err_wrps_lin  / va,
        err_b_poly    / va,
        err_b_lin     / va
      ))
    }
    
    
    # ============================================================
    # CALCOLO ERRORI
    # ============================================================
    
    # ------------------------------------------------------------
    # DGP POLINOMIALE
    # ------------------------------------------------------------
    
    errs_poly <- calc_errors(
      A_poly = A1_poly,
      A_lin  = A1_lin,
      y_true = y1,
      x_poly = x_p,
      x_lin  = x_l
    )
    
    
    # ------------------------------------------------------------
    # DGP LINEARE
    # ------------------------------------------------------------
    
    errs_lin <- calc_errors(
      A_poly = A2_poly,
      A_lin  = A2_lin,
      y_true = y2,
      x_poly = x_p,
      x_lin  = x_l
    )
    
    
    # ============================================================
    # SALVATAGGIO RISULTATI
    # ============================================================
    
    result[it, 1:4] <- errs_poly
    result[it, 5:8] <- errs_lin
    
    
    # ============================================================
    # WINS
    # ============================================================
    
    # ------------------------------------------------------------
    # POLY vs LINEAR sul DGP POLINOMIALE
    # ------------------------------------------------------------
    
    result[it, 9:10] <- as.numeric(
      errs_poly[c(1,3)] < errs_poly[c(2,4)]
    )
    
    
    # ------------------------------------------------------------
    # POLY vs LINEAR sul DGP LINEARE
    # ------------------------------------------------------------
    
    result[it, 11:12] <- as.numeric(
      errs_lin[c(1,3)] < errs_lin[c(2,4)]
    )
    
    
    # ------------------------------------------------------------
    # TRUE DGP
    # Come nella vecchia simulazione:
    # POLY sul DGP POLY vs LIN sul DGP LIN
    # ------------------------------------------------------------
    
    result[it, 13:14] <- as.numeric(
      errs_poly[c(1,3)] < errs_lin[c(2,4)]
    )
    
    # ============================================================
    # FINE ITERAZIONE
    # ============================================================
    
  }
  
  
  # ============================================================
  # RISULTATI FINALI
  # ============================================================
  
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
    
    win_rates = colMeans(
      result[, 9:14],
      na.rm = TRUE
    ),
    
    raw_results = result
  )
}    

s2 <- run_simulation2(
  iter=500,
  N = 100,
  Cx0 = 5,
  Cy0 = 7,
  response_type="ordinal",
  distance="wd"
)
s2$win_rates  

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

run_simulation_caseI2 <- function(
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
  
  result <- matrix(NA, nrow = iter, ncol = 12)
  
  colnames(result) <- c(
    "W1_LIN",
    "SCS_LIN",
    "CVM_LIN",
    "R2_W1_LIN",
    "R2_SCS_LIN",
    "R2_CVM_LIN",
    
    "W1_POLY",
    "SCS_POLY",
    "CVM_POLY",
    "R2_W1_POLY",
    "R2_SCS_POLY",
    "R2_CVM_POLY"
  )
  
  
  for(it in 1:iter){
    
    # =====================================================
    # GENERAZIONE DGP
    # =====================================================
    
    data <- generate_linear_poly_data(
      N = N,
      Cx = Cx0,
      Cy = Cy0,
      degree = degree,
      response_type = response_type,
      lambda = lambda,
      alpha_dir = alpha_dir
    )
    
    
    # =====================================================
    # DGP LINEARE
    # =====================================================
    
    X_lin <- data$X
    Y_lin <- data$Y_linear
    
    
    # =====================================================
    # DGP POLINOMIALE
    # =====================================================
    
    X_poly <- data$Z
    Y_poly <- data$Y_poly
    
    
    # =====================================================
    # WEIGHTS
    # =====================================================
    
    weights_unit <- rep(1, Cy0 - 1)
    
    wsum_u <- sum(weights_unit)
    wsum_cvm <- sum(weights_unit^2)
    
    
    # =====================================================
    # =====================================================
    #                DGP LINEARE
    # =====================================================
    # =====================================================
    
    
    # =====================================================
    # LISTE
    # =====================================================
    
    X_list_lin <- lapply(
      1:N,
      function(i) X_lin[i,]
    )
    
    Y_list_lin <- lapply(
      1:N,
      function(i) Y_lin[i,]
    )
    
    
    # =====================================================
    # W1 - LINEARE
    # =====================================================
    
    A_unit_lin <- solve_simplex_lp(
      X_list_lin,
      Y_list_lin,
      weights_unit
    )$A
    
    
    pred_u_lin <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_u_lin[i,] <- as.vector(
        A_unit_lin %*% X_lin[i,]
      )
    }
    
    
    SSE_u_lin <- 0
    
    for(i in 1:N){
      SSE_u_lin <- SSE_u_lin +
        wd(
          weights_unit,
          Y_lin[i,],
          pred_u_lin[i,]
        )
    }
    
    result[it,"W1_LIN"] <-
      SSE_u_lin / (wsum_u * N)
    
    
    result[it,"R2_W1_LIN"] <- compute_R2(
      Y_lin,
      X_lin,
      A_unit_lin,
      distance = "wd",
      weights = weights_unit
    )$R2
    
    
    # =====================================================
    # SCS - LINEARE
    # =====================================================
    
    A_scs_lin <- fit_simplex_nominal(
      X_lin,
      Y_lin,
      degree = 1
    )$A
    
    
    pred_scs_lin <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_scs_lin[i,] <- as.vector(
        A_scs_lin %*% X_lin[i,]
      )
    }
    
    
    SSE_scs_lin <- 0
    
    for(i in 1:N){
      SSE_scs_lin <- SSE_scs_lin +
        sym_chi2(
          Y_lin[i,],
          pred_scs_lin[i,]
        )
    }
    
    result[it,"SCS_LIN"] <-
      SSE_scs_lin / (4 * N)
    
    
    result[it,"R2_SCS_LIN"] <- compute_R2(
      Y_lin,
      X_lin,
      A_scs_lin,
      distance = "scs"
    )$R2
    
    
    # =====================================================
    # CVM - LINEARE
    # =====================================================
    
    A_cvm_lin <- fit_simplex_ordinal(
      X_lin,
      Y_lin,
      degree = 1,
      distance = "cvm"
    )$A
    
    
    pred_cvm_lin <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_cvm_lin[i,] <- as.vector(
        A_cvm_lin %*% X_lin[i,]
      )
    }
    
    
    SSE_cvm_lin <- 0
    
    for(i in 1:N){
      SSE_cvm_lin <- SSE_cvm_lin +
        cvm_distance(
          weights_unit,
          Y_lin[i,],
          pred_cvm_lin[i,]
        )
    }
    
    result[it,"CVM_LIN"] <-
      SSE_cvm_lin / (wsum_cvm * N)
    
    
    result[it,"R2_CVM_LIN"] <- compute_R2(
      Y_lin,
      X_lin,
      A_cvm_lin,
      distance = "cvm",
      weights = weights_unit
    )$R2
    
    
    # =====================================================
    # =====================================================
    #                DGP POLINOMIALE
    # =====================================================
    # =====================================================
    
    
    # =====================================================
    # LISTE
    # =====================================================
    
    X_list_poly <- lapply(
      1:N,
      function(i) X_poly[i,]
    )
    
    Y_list_poly <- lapply(
      1:N,
      function(i) Y_poly[i,]
    )
    
    
    # =====================================================
    # W1 - POLINOMIALE
    # =====================================================
    
    A_unit_poly <- solve_simplex_lp(
      X_list_poly,
      Y_list_poly,
      weights_unit
    )$A
    
    
    pred_u_poly <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_u_poly[i,] <- as.vector(
        A_unit_poly %*% X_poly[i,]
      )
    }
    
    
    SSE_u_poly <- 0
    
    for(i in 1:N){
      SSE_u_poly <- SSE_u_poly +
        wd(
          weights_unit,
          Y_poly[i,],
          pred_u_poly[i,]
        )
    }
    
    result[it,"W1_POLY"] <-
      SSE_u_poly / (wsum_u * N)
    
    
    result[it,"R2_W1_POLY"] <- compute_R2(
      Y_poly,
      X_poly,
      A_unit_poly,
      distance = "wd",
      weights = weights_unit
    )$R2
    
    
    # =====================================================
    # SCS - POLINOMIALE
    # =====================================================
    
    A_scs_poly <- fit_simplex_nominal(
      X_poly,
      Y_poly,
      degree = 1
    )$A
    
    
    pred_scs_poly <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_scs_poly[i,] <- as.vector(
        A_scs_poly %*% X_poly[i,]
      )
    }
    
    
    SSE_scs_poly <- 0
    
    for(i in 1:N){
      SSE_scs_poly <- SSE_scs_poly +
        sym_chi2(
          Y_poly[i,],
          pred_scs_poly[i,]
        )
    }
    
    result[it,"SCS_POLY"] <-
      SSE_scs_poly / (4 * N)
    
    
    result[it,"R2_SCS_POLY"] <- compute_R2(
      Y_poly,
      X_poly,
      A_scs_poly,
      distance = "scs"
    )$R2
    
    
    # =====================================================
    # CVM - POLINOMIALE
    # =====================================================
    
    A_cvm_poly <- fit_simplex_ordinal(
      X_poly,
      Y_poly,
      degree = 1,
      distance = "cvm"
    )$A
    
    
    pred_cvm_poly <- matrix(
      0,
      nrow = N,
      ncol = Cy0
    )
    
    for(i in 1:N){
      pred_cvm_poly[i,] <- as.vector(
        A_cvm_poly %*% X_poly[i,]
      )
    }
    
    
    SSE_cvm_poly <- 0
    
    for(i in 1:N){
      SSE_cvm_poly <- SSE_cvm_poly +
        cvm_distance(
          weights_unit,
          Y_poly[i,],
          pred_cvm_poly[i,]
        )
    }
    
    result[it,"CVM_POLY"] <-
      SSE_cvm_poly / (wsum_cvm * N)
    
    
    result[it,"R2_CVM_POLY"] <- compute_R2(
      Y_poly,
      X_poly,
      A_cvm_poly,
      distance = "cvm",
      weights = weights_unit
    )$R2
    
    
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
sI <- run_simulation_caseI2(
  iter=500,
  N = 100,
  Cx = 5,
  Cy = 5,
  response_type = "ordinal"
)
