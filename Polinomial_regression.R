#install.packages("CVXR",type = "binary")

library(CVXR)

# Recursive function to generate multi-indices that sum to degree k
generate_multiindices <-function(n_vars, degree) {
 if (n_vars == 1) return(matrix(degree, nrow = 1))
  res <-matrix(nrow = 0, ncol = n_vars)
  for (i in 0:degree) {
   sub_res <-generate_multiindices(n_vars- 1, degree- i)
   res <-rbind(res, cbind(i, sub_res))
  }
 return(res)
}


#generate_multiindices(3,2)

 # Function to map predictors X into the Bernstein basis Z^(k)
bernstein_basis <-function(X, degree, alpha_weights=NULL) {
   N <-nrow(X)
   n_plus_1 <-ncol(X)
   n <- n_plus_1 - 1
   indices <-generate_multiindices(n_plus_1, degree)
   D <-nrow(indices)
  
   Z <-matrix(0, nrow = N, ncol = D)
   for (d in 1:D) {
     gamma <-indices[d, ]
     # Multinomial coefficient: k! / (gamma_1! ... gamma_{n+1}!)
     coef <-factorial(degree) / prod(factorial(gamma))
     Z[, d] <-coef * apply(X, 1, function(x) prod(x^gamma))
   }
   # Default to equidistant predictor vertices ( alpha _ i = 1) -> handles nominal   data gracefully
    if ( is.null ( alpha_weights ) ) {
      alpha_weights <- rep (1 , n )
      }
   
    if ( length ( alpha_weights ) != n || any ( alpha_weights <= 0) ) {
      stop ( sprintf ( " weights must be a strictly positive vector of length % d. " ,n ) )
      }
   
    # Compute cumulative ground costs : c _ 1 = 0 , c _ i = sum ( alpha _ {1:( i -1) })
    c_costs <- c (0 , cumsum ( alpha_weights ) )
   
    # Compute the polynomial cost C ( gamma )
    gamma_costs <- apply ( indices , 1 , function ( gamma ) sum ( gamma * c_costs ) )
   
    # Sort indices by cost ( the ’ order ’ function is stable , preservinglexicographical order for ties )
 sort_idx <- order ( gamma_costs )
 Z <- Z [ , sort_idx , drop = FALSE ]
   return(Z)
}

#bernstein_basis(x,2)

fit_simplex_ordinal <-function(X, Y, degree = 2, weights = NULL, alpha_weights=NULL) {
   N <-nrow(Y)
   m_plus_1 <-ncol(Y)
   m <-m_plus_1- 1
   Z <-bernstein_basis(X, degree,alpha_weights)
   D <-ncol(Z)
  
   # Default to equidistant vertices (a_j = 1) if weights are not provided
   if (is.null(weights)) {
     a_weights <-rep(1, m)
     }
   if (length(weights) != m || any(weights <= 0)) {
     stop(sprintf("a_weights must be a strictly positive vector of length %d.", m))
     }
  
   # Truncated accumulation matrix L_m
   Lm <-matrix(0, nrow = m, ncol = m_plus_1)
   for (i in 1:m) {
     Lm[i, 1:i] <-1
     }
  
   # Diagonal weight matrix for ground distances
   W_diag <-diag(weights, nrow = m, ncol = m)
  
   A_var <-Variable(m_plus_1, D)
   Y_pred <-Z %*% t(A_var)
  
   # Map PMFs to CDFs and apply the geometric ground distances
   # (Post-multiplying by W_diag scales each column j by a_j)
   weighted_F_true <-(Y %*% t(Lm)) %*% W_diag
   weighted_F_pred <-(Y_pred %*% t(Lm)) %*% W_diag
  
   # Single-step strictly convex optimization
   # CVXR’s sum_squares will automatically square the a_j coefficient yielding a_j^2
   objective <-Minimize(sum_squares(weighted_F_true-weighted_F_pred))
  
   constraints <-list(
     A_var >= 0,
     matrix(1, 1, m_plus_1) %*% A_var == matrix(1, 1, D)
     )
  
   prob <-Problem(objective, constraints)
   res <-solve(prob, solver = "ECOS")
  
   A_opt <-res$getValue(A_var)
  
   # Clean up numerical rounding errors
   A_opt[A_opt < 0] <-0
   A_opt <-sweep(A_opt, 2, colSums(A_opt), "/")
   return(list(A = A_opt, fitted = Z %*% t(A_opt)))
}   

res<-fit_simplex_ordinal(x,y)
res$A
RMSE <- sqrt(mean((y - res$fitted)^2))
RMSE
RMSE_obs <- sqrt(rowMeans((y - res$fitted)^2))
RMSE_obs

##############con regolarizzazione##################

# ==========================================================
# Genera gli archi della griglia Bernstein
# ==========================================================

generate_edges <- function(indices){
  
  D <- nrow(indices)
  
  edges <- list()
  counter <- 1
    
    if(D > 1){
      for(i in 1:(D-1)){
        edges[[counter]] <- c(i, i+1)
        counter <- counter + 1
      }
    }
    
    return(edges)
  
}



# ==========================================================
# Genera triple per differenze seconde
# ==========================================================

generate_second_order <- function(indices){
  
  D <- nrow(indices)
  p <- ncol(indices)
  
  triples <- list()
  counter <- 1
  
    
    if(D > 2){
      
      for(i in 1:(D-2)){
        
        triples[[counter]] <- c(i, i+1, i+2)
        counter <- counter + 1
        
      }
      
    }
    
    return(triples)
    
}




# ==========================================================
# Modello simplex ordinale Bernstein con regolarizzazione
# ==========================================================


fit_simplex_ordinal<-function(
  X,
  Y,
  degree = 2,
  distance = c("cvm","wasserstein"),
    weights = NULL,
  alpha_weights=NULL,
    lambda1 = 0,
    lambda2 = 0
){
  
  
  N <- nrow(Y)
  
  m_plus_1 <- ncol(Y)
  
  m <- m_plus_1 - 1
  
  
  
  # --------------------------
  # Bernstein basis
  # --------------------------
  
  Z <- bernstein_basis(
    X,
    degree, alpha_weights
  )
  
  indices <- generate_multiindices(
    ncol(X),
    degree
  )
  
  
  D <- ncol(Z)
  
  # --------------------------
  # weights
  # --------------------------
  
  if(is.null(weights)){
    weights <- rep(1,m)
  }
  
  
  if(length(weights)!=m ||
     any(weights<=0)){
    
    stop(
      "a_weights must be positive"
    )
  }
  
  
  
  # --------------------------
  # CDF transformation
  # --------------------------
  
  Lm <- matrix(
    0,
    nrow=m,
    ncol=m_plus_1
  )
  
  
  for(i in 1:m){
    Lm[i,1:i] <- 1
  }
  
  
  W <- diag(weights)
  
  
  # --------------------------
  # Variable
  # --------------------------
  
  A <- Variable(
    m_plus_1,
    D
  )
  
  
  Y_hat <- Z %*% t(A)
  
  F_true <- 
    (Y %*% t(Lm)) %*% W
  
  F_pred <-
    (Y_hat %*% t(Lm)) %*% W
  
  
  # --------------------------
  # Loss
  # --------------------------
  
  distance <- match.arg(distance)
  
  if(distance=="cvm"){
    
    loss <- sum_squares(F_true - F_pred)
    
  }else{
    
    #weights <- matrix(  a_weights,nrow=N,ncol=m,  byrow=TRUE)
    #loss <- sum(multiply(  weights,abs(F_true - F_pred)))
    loss <- sum(
      abs(F_true - F_pred)
    )
  }
  
  
  # --------------------------
  # Regolarizzazione primo ordine
  # --------------------------
  
  reg1 <- 0
  
  
  if(lambda1>0){
    
    edges <- generate_edges(
      indices
    )
    
    
    for(e in edges){
      
      reg1 <- reg1 +
        sum_squares(
          A[,e[1]]-
            A[,e[2]]
        )
    }
    
  }
  
  
  # --------------------------
  # Regolarizzazione secondo ordine
  # --------------------------
  
  reg2 <- 0
  
  
  if(lambda2>0){
    
    triples <- generate_second_order(
      indices
    )
    
    
    for(t in triples){
      
      reg2 <- reg2 +
        sum_squares(
          A[,t[1]]
          -2*A[,t[2]]
          +A[,t[3]]
        )
    }
    
  }
  
  
  # --------------------------
  # Objective
  # --------------------------
  
  objective <- Minimize(
    loss +
      lambda1*reg1 +
      lambda2*reg2
  )
  
  
  # --------------------------
  # Constraints
  # --------------------------
  
  constraints <- list(
    
    A >= 0,
    
    matrix(
      1,
      1,
      m_plus_1
    ) %*% A
    ==
      matrix(
        1,
        1,
        D
      )
    
  )
  
  
  # --------------------------
  # Solve
  # --------------------------
  
  problem <- Problem(
    objective,
    constraints
  )
  
  
  result <- solve(
    problem,
    solver="ECOS"
  )
  
  
  A_opt <- result$getValue(A)
  
  
  # --------------------------
  # Pulizia numerica
  # --------------------------
  
  A_opt[A_opt<0] <- 0
  
  
  A_opt <- sweep(  A_opt,  2,  colSums(A_opt),  "/")
  
  
  return(
    list(
      A=A_opt,
      fitted=Z %*% t(A_opt),
      objective=result$value,
      indices=indices
    )
  )
  
}


select_lambda_new <- function(
    X,
    Y,
    degree = 2,
    weights = NULL,
    alpha_weights=NULL,
    lambda1_grid,
    lambda2_grid,
    method = c("cv", "gcv"),
    distance = c("cvm","wasserstein"),
    K = NULL,
    seed = 123,
    verbose = FALSE
){
  distance <- match.arg(distance)
  method <- match.arg(method)
  
  
  N <- nrow(X)
  
  m <- ncol(Y)-1
  
  if(is.null(weights)){
    weights <- rep(1,m)
  }
  
  if(length(weights)!=m || any(weights<=0)){
    stop(
      "a_weights must be a positive vector of length ncol(Y)-1"
    )
  }
  
  lambda1_grid <- sort(unique(lambda1_grid))
  lambda2_grid <- sort(unique(lambda2_grid))
  
  
  # =====================================================
  # Funzione errore quadratico sulle CDF
  # =====================================================
  
  compute_fit_single <- function(A, X_new, Y_new){
    
    Z_new <- bernstein_basis(
      X_new,
      degree,alpha_weights
    )
    
    Y_hat <- Z_new %*% t(A)
    
    m_plus_1 <- ncol(Y_new)
    m <- m_plus_1 - 1
    
    Lm <- matrix(
      0,
      nrow = m,
      ncol = m_plus_1
    )
    
    for(i in 1:m){
      Lm[i,1:i] <- 1
    }
    
    F_true <- Y_new %*% t(Lm)
    F_pred <- Y_hat %*% t(Lm)
    
    if(distance == "cvm"){
      
      W <- diag(weights)
      
      err <- sum(
        (((F_true - F_pred) %*% W)^2)
      )
      
    } else {
      
      err <- sum(
        abs(F_true - F_pred) *
          matrix(
            weights,
            nrow = nrow(F_true),
            ncol = m,
            byrow = TRUE
          )
      )
      
    }
    
    if(is.na(err) || is.infinite(err)){
      return(Inf)
    }
    
    err
  }
  
  
  compute_fit <- function(A, X_all, Y_all){
    
    Z_all <- bernstein_basis(
      X_all,
      degree,alpha_weights
    )
    
    Y_hat <- Z_all %*% t(A)
    
    m_plus_1 <- ncol(Y_all)
    m <- m_plus_1 - 1
    
    Lm <- matrix(
      0,
      nrow = m,
      ncol = m_plus_1
    )
    
    for(i in 1:m){
      Lm[i,1:i] <- 1
    }
    
    F_true <- Y_all %*% t(Lm)
    F_pred <- Y_hat %*% t(Lm)
    
    if(distance == "cvm"){
      
      W <- diag(weights)
      
      fit <- sum(
        (((F_true - F_pred) %*% W)^2)
      )
      
    } else {
      
      fit <- sum(
        abs(F_true - F_pred) *
          matrix(
            weights,
            nrow = nrow(F_true),
            ncol = m,
            byrow = TRUE
          )
      )
      
    }
    
    fit
  }
  
  
  
  # =====================================================
  # CROSS VALIDATION
  # =====================================================
  
  if(method=="cv"){
    
    
    if(is.null(K)){
      K <- N
    }
    
    
    if(K>N){
      stop("K must be lower than number of observations")
    }
    
    
    set.seed(seed)
    
    
    if(K==N){
      
      folds <- seq_len(N)
      
    } else {
      
      folds <- sample(
        rep(1:K,length.out=N)
      )
      
    }
    
    
    
    cv_values <- matrix(
      Inf,
      length(lambda1_grid),
      length(lambda2_grid)
    )
    
    
    rownames(cv_values) <-
      paste0("lam1_",lambda1_grid)
    
    colnames(cv_values) <-
      paste0("lam2_",lambda2_grid)
    
    
    for(i1 in seq_along(lambda1_grid)){
      
      
      lam1 <- lambda1_grid[i1]
      
      
      for(i2 in seq_along(lambda2_grid)){
        
        
        lam2 <- lambda2_grid[i2]
        
        
        if(verbose){
          cat(
            "\nlambda1 =",lam1,
            " lambda2 =",lam2,"\n"
          )
        }
        
        
        cv_error <- 0
        
        valid <- TRUE
        
        
        for(k in unique(folds)){
          
          
          train_idx <- which(folds!=k)
          
          test_idx <- which(folds==k)
          
          
          res <- tryCatch({
            
            fit_simplex_ordinal(
              X[train_idx,,drop=FALSE],
              Y[train_idx,,drop=FALSE],
              degree = degree,
              distance = distance,
              weights = weights,
              lambda1 = lam1,
              lambda2 = lam2
            )
            
          },
          error = function(e) NULL)
          
          
          if(is.null(res)){
            
            valid <- FALSE
            break
            
          }
          
          
          A_hat <- res$A
          
          
          for(t in test_idx){
            
            
            err <- compute_fit_single(
              A_hat,
              X[t,,drop=FALSE],
              Y[t,,drop=FALSE]
            )
            
            
            if(is.infinite(err)){
              
              valid <- FALSE
              break
              
            }
            
            
            cv_error <- cv_error + err
            
            
          }
          
          
          if(!valid) break
          
          
        }
        
        
        
        if(valid){
          
          cv_values[i1,i2] <- cv_error
          
          
          if(verbose)
            cat("CV error =",cv_error,"\n")
          
          
        }
        
      }
      
    }
    
    
    
    if(all(is.infinite(cv_values))){
      stop("All configurations discarded")
    }
    
    
    
    idx <- which(
      cv_values==min(cv_values),
      arr.ind=TRUE
    )[1,]
    
    
    return(list(
      
      method="cv",
      
      best_lambda1=lambda1_grid[idx[1]],
      
      best_lambda2=lambda2_grid[idx[2]],
      
      cv_values=cv_values,
      
      lambda1_grid=lambda1_grid,
      
      lambda2_grid=lambda2_grid,
      
      K=K,
      
      seed=seed
      
    ))
    
    
  }
  
  
  # =====================================================
  # GCV
  # =====================================================
  
  if(method=="gcv"){
    
    
    n1 <- length(lambda1_grid)
    
    n2 <- length(lambda2_grid)
    
    
    fit_values <- rep(Inf,n1*n2)
    
    score_values <- rep(Inf,n1*n2)
    
    
    cc <- 0
    
    
    
    for(i1 in seq_along(lambda1_grid)){
      
      
      for(i2 in seq_along(lambda2_grid)){
        
        
        cc <- cc+1
        
        
        lam1 <- lambda1_grid[i1]
        
        lam2 <- lambda2_grid[i2]
        
        
        
        res <- tryCatch(
          
          
          fit_simplex_ordinal(
            X,
            Y,
            degree = degree,
            distance = distance,
            weights = weights,
            lambda1 = lam1,
            lambda2 = lam2
          ),
          
          
          error=function(e) NULL)
        
        
        
        if(is.null(res))
          next
        
        
        
        A_hat <- res$A
        
        
        
        fit_values[cc] <-
          compute_fit(
            A_hat,
            X,
            Y
          )
        
        
        
        # stessa idea del tuo GCV
        
        sv <- svd(A_hat)$d
        
        sv <- sv[sv>1e-6]
        
        
        df <- sum(
          sv^2/(sv^2+lam1+lam2)
        )
        
        
        
        n <- N*(ncol(Y)-1)
        
        
        score_values[cc] <-
          fit_values[cc]/
          (max(n-df,1)^2)
        
        
        
      }
      
    }
    
    
    
    best <- which.min(score_values)
    
    
    
    i1 <- (best-1)%/%n2+1
    
    i2 <- (best-1)%%n2+1
    
    
    
    return(list(
      
      method="gcv",
      
      best_lambda1=lambda1_grid[i1],
      
      best_lambda2=lambda2_grid[i2],
      
      score_values=score_values,
      
      fit_values=fit_values,
      
      lambda1_grid=lambda1_grid,
      
      lambda2_grid=lambda2_grid
      
    ))
    
    
  }
  
  
}

fit_simplex_nominal <- function (X , Y , degree = 2 , alpha_weights = NULL , max_iter =50 , tol = 1e-5) {
   N <- nrow ( Y )
   m_plus_1 <- ncol ( Y )

   # Map and topologically sort the predictors using alpha_weights
   Z <- bernstein_basis (X , degree ,alpha_weights )
   D <- ncol ( Z )
   
    # Initialization : constant prediction equal to the marginal mean
    Y_hat <- matrix ( colMeans ( Y ) , nrow = N , ncol = m_plus_1 , byrow = TRUE )
    A_prev <- matrix (1/m_plus_1 , nrow = m_plus_1 , ncol = D )
   
    for ( iter in 1: max_iter ) {
      # Numerical clipping : Max weight restricted to 10^6 to prevent Hessian overflow
      denom <- pmax (( Y + Y_hat ) / 2 , 1e-6)
      W <- 1 / denom
     
      A_var <- Variable ( m_plus_1 , D )
      Y_pred <- Z %*% t ( A_var )
     
      # Weighted least squares inner loop
      objective <- Minimize ( sum_squares ( sqrt( W ) * ( Y - Y_pred ) ) )
     
      # Strict simplex constraints for the parameter matrix A
      constraints <- list (
        A_var >= 0 ,
        matrix(1 , 1 , m_plus_1) %*% A_var == matrix (1 , 1 , D )
        )
     
      prob <- Problem( objective , constraints )
      res <- solve( prob , solver = "ECOS" )
     
      A_new <- res$getValue( A_var )
      Y_hat <- Z %*% t( A_new )
     
      if ( max ( abs( A_new - A_prev ) ) < tol ) break
      A_prev <- A_new
      }
   
    # Clean up final numerical rounding errors
    A_new [ A_new < 0] <- 0
    A_new <- sweep ( A_new , 2 , colSums ( A_new ) , "/" )
   
    return ( list ( A = A_new , fitted = Y_hat ) )
  }

  
library(codalm)
data("educFM")
father <- as.matrix(educFM[,2:4])
ya <- father/rowSums(father)
mother <- as.matrix(educFM[,5:7])
xa <- mother/rowSums(mother)
x<-matrix(0,nrow=dim(ya)[1],ncol=dim(ya)[2])
x[,1]<-xa[,3]
x[,2]<-xa[,2]
x[,3]<-xa[,1]
y<-matrix(0,nrow=dim(ya)[1],ncol=dim(ya)[2])
y[,1]<-ya[,3]
y[,2]<-ya[,2]
y[,3]<-ya[,1]
weights<-c(1,1)

lambda_grid <- c(0, 0.005, 0.01, 0.02, 0.05, 0.08, 0.1)
library(OrdinalCompositions)
N<-length(x[,1])
ydata <- split(y, seq_len(N))
xdata <- split(x, seq_len(N))
xdatabernstein<-split(bernstein_basis(x,2), seq_len(N))

res<-select_lambda_new(x,y,degree=2,lambda1_grid=lambda_grid,lambda2_grid=lambda_grid,method="gcv",distance="cvm")
fit<-fit_simplex_ordinal(x,y,degree=2,distance="cvm",lambda1=res$best_lambda1,lambda2=res$best_lambda2)
fit$A
compute_R2(y, bernstein_basis(x,2), fit$A, weights)

res<-select_lambda(xdatabernstein,ydata,weights,lambda1_grid=lambda_grid,lambda2_grid=lambda_grid,method="cv")
fit<-solve_simplex_lp(xdatabernstein,ydata,weights,lambda1=res$best_lambda1,lambda2=res$best_lambda2)
fit$A
compute_R2(y, bernstein_basis(x,2), fit$A, weights)

t(codalm(y,bernstein_basis(x,2)))

fit<-fit_simplex_nominal(x,y,degree=2,weights)
fit$A
compute_R2(y, bernstein_basis(x,2), fit$A, weights)


ordinal_regression_simplex(
  split(bernstein_basis(x,2), seq_len(N)),
  ydata,
  weights = weights,
  lambda1_grid = 0,
  lambda2_grid = 0,
  method = "gcv",
  do_bootstrap = FALSE,
  B = 200,
  compute_opi = TRUE,
  compute_R2 = TRUE,
  compute_OCC = TRUE
)
