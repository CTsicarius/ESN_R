generate_W <- function(Nx, density, sigma = 1, ro = 1){
  #Matrix
  W <- Matrix::rsparsematrix(Nx, Nx, density, rand.x = rnorm)*sigma
  vaps <- eigen(W, only.values = TRUE)$values
  W <- ro*W/max(abs(vaps))
  return(W)
}

generate_Win <- function(Nx, Nu, mu = 0, sigma = 1){
  Win <- matrix(rnorm(Nx * (1 + Nu), mu, sigma), Nx, 1 + Nu)
  return(Win)
}

calculate_x <- function(u, W, Win, alpha){
  Nx <- dim(Win)[1]
  T0 <- dim(u)[2]
  x <- matrix(0, Nx, T0)
  #x_bar <- as.matrix(tanh(Win %*% t(cbind(1, t(u[, 1])))))
  x_bar <- tanh(Win %*% t(cbind(1, t(u[, 1]))))
  x[, 1] <- alpha * x_bar
  if(T0 > 1){
    for(i in 2:T0){
      x_bar <- as.matrix(tanh(Win %*% t(cbind(1, t(u[, i]))) + W %*% x[, i - 1]))
      x[, i] <- (1 - alpha) * x[, i - 1] + alpha * x_bar
    }
  }
  return(x)  
}

cppFunction("arma::mat calculate_x_c(arma::mat u, arma::mat W, arma::mat Win, double alpha) {
                        int Nx = Win.n_rows;
                        int T0 = u.n_cols;
                        arma::mat x = arma::zeros(Nx, T0);
                        arma::mat x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(0)));
                        //HASTA AQUI TODO BIEN
                        x.col(0) = alpha*x_bar;
                        for(int i = 1; i < T0; ++i){
                            x_bar = arma::tanh(Win*arma::join_cols(arma::ones(1, 1), u.col(i)) + W*x.col(i - 1));
                            x.col(i) = (1. - alpha)*x.col(i - 1) + alpha*x_bar;
                        }
                        return(x);
            }", depends='RcppArmadillo')

train_Wout <- function(y, u, W, Win, x = NULL, alpha, beta){
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  Ny <- dim(y)[1]
  Nu <- dim(u)[1]
  Nx <- dim(x)[1]
  T0 <- dim(u)[2]
  input <- rbind(matrix(1, 1, T0), rbind(u, x))
  Wout <- t(solve(input %*% t(input) + beta*diag(Nx + Nu + 1), input %*% t(y)))
  return(Wout)
}


solve_train_Wout_all <- function(y, u, W, Win, alpha){
  #CVXR
  x <- calculate_x(u, W, Win, alpha)
  Ny <- dim(y)[1]
  Nu <- dim(u)[1]
  Nx <- dim(x)[1]
  T0 <- dim(u)[2]
  input <- rbind(matrix(1, 1, T0), rbind(u, x))
  Wout <- CVXR::Variable(Ny, 1 + Nu + Nx)
  objective <- CVXR::Minimize(sum((y - (Wout %*% input))^2) / T0)
  problem <- CVXR::Problem(objective)
  return(solve(problem))
}

make_esn_predictions <- function(u, W, Win, x = NULL, alpha, Wout){
  if(is.null(x)){
    x <- calculate_x(u, W, Win, alpha)
  } 
  T0 <- dim(u)[2]
  return(Wout %*% rbind(matrix(1, 1, T0), rbind(u, x)))
}





