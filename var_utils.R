generate_var1 <- function(T0, phi0, PHI, sigma){
  #mvrnorm package MASS
  N <- dim(PHI)[1]
  r_mat <- matrix(0, nrow = N, ncol = T0)
  r_mat[,1] <- rep(1, N)
  for (i in 2:T0) {
    r_mat[, i] <- phi0 + PHI%*%r_mat[, i-1] + MASS::mvrnorm(1, rep(0, N), sigma)
  }
  return(r_mat)
}

generate_PHI_normal <- function(N, mu = 0, sigma = 1, ro = 1){
  PHI <- matrix(rnorm(N*N, mu, sigma), N, N)
  vaps <- eigen(PHI, only.values = TRUE)$values
  PHI <- ro*PHI/max(abs(vaps))
  return(PHI)
}

estimate_parameters <- function(data){
  #CVXR
  N <- dim(data)[1]
  T0 <- dim(data)[2]
  data1 <- data[,2:T0]
  data0 <- data[,1:T0 - 1]
  data0 <- rbind(matrix(1, 1, T0 - 1), data0)
  PHI_hat <- t(solve(data0 %*% t(data0), data0 %*% t(data1)))
  return(PHI_hat)
}

estimate_all <- function(data){
  #CVXR
  N <- dim(data)[1]
  T0 <- dim(data)[2]
  data1 <- data[,2:T0]
  data0 <- data[,1:T0 - 1]
  data0 <- rbind(matrix(1, 1, T0 - 1), data0)
  PHI_hat <- CVXR::Variable(N, N + 1, name = 'PHIII')
  objective <- CVXR::Minimize( sum( (data1 - (PHI_hat%*%data0))^2) /(T0 - 1))
  problem <- CVXR::Problem(objective)
  return(solve(problem))
}

make_var_predictions <- function(data, PHI){
  N <- dim(data)[1]
  T0 <- dim(data)[2]
  data1 <- data[,2:T0]
  data0 <- data[,1:T0 - 1]
  data0 <- rbind(matrix(1, 1, T0 - 1), data0)
  res <-  PHI%*%data0
  return(res)
}

MSE_error <- function(real_values, pred_values){
  T0 <- dim(real_values)[2]
  return(sum((real_values - pred_values)^2)/T0)
}
