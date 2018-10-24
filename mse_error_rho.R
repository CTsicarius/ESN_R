library(MASS)
library(Rcpp)
library(RcppArmadillo)
wd <- getwd()
source(paste(wd, 'var_utils.R', sep = '/'))
source(paste(wd, 'esn_utils.R', sep = '/'))
#SIM PARAMETERS
N_sim <- 1
print_det = FALSE
#DATA PARAMETERS
seed = 42
N <- 1
Ny <- 1
Nu <- 1
PHI_ro = 0.2
N_test <- 10000
sigma <- 0.1
#ESN PARAMETERS
Nx <- 5000
w_ro <- 0.001
w_dens <- 0.5
alpha <- 1
beta <- 5

set.seed(42)
number_of_train <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000)
#rho_array <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
rho_array <- c(0, 0.25, 0.5, 0.75, 0.9)

MSE_array_var <- matrix(0, length(number_of_train), 1)
MSE_array_esn <- matrix(0, length(number_of_train), length(rho_array))
phi <- matrix(rnorm(1:N), N, 1)
PHI1 <- generate_PHI_normal(N = N, mu = 0, sigma = 1, ro = 0.8)
PHI2 <- generate_PHI_normal(N = N, mu = 0, sigma = 1, ro = 0.95)
#global_train_data <-generate_var2(number_of_train[length(number_of_train)], phi, PHI1, PHI2, diag(sigma, N, N))#4.5 secs
#global_train_data <- generate_var1(number_of_train[length(number_of_train)], phi, PHI1, diag(sigma, N, N)) #4.5 secs
global_train_data = matrix(c(1, 1, 1, 0, 0, 1, 1, 0), 1, 50000)
#test_data <- generate_var2(N_test, phi, PHI1, PHI2, diag(sigma, N, N))
#test_data <- generate_var1(N_test, phi, PHI1, diag(sigma, N, N))
test_data =  matrix(c(1, 0, 0, 1, 1, 0, 1, 1), 1, 10000)
time0 <- proc.time()
print('OK1')
Win <- generate_Win(Nx = Nx, Nu = Nu, mu = 0, sigma = 1)
for (rho_idx in 1:length(rho_array)){
  W <- as.matrix(generate_W(Nx = Nx, density = w_dens, sigma = 1, ro = rho_array[rho_idx]))
  train_x <- calculate_x_c(u = matrix(global_train_data[, 1:number_of_train[length(number_of_train)]], 1, number_of_train[length(number_of_train)]), 
                           W = W, Win = Win, alpha = alpha)
  test_x <- calculate_x_c(u = test_data, W = W, Win = Win, alpha = alpha)
  for (i in 1:length(number_of_train)){
    #CODE FOR VAR ERROR
    train_data <- matrix(global_train_data[, 1:number_of_train[i]], 1, number_of_train[i])
    PHI_est <- estimate_var2_parameters(train_data)
    var_pred_val <- make_var2_predictions(test_data, PHI_est)
    MSE_array_var[i,  1] <- MSE_error(test_data[, 3:N_test, drop = FALSE], var_pred_val)
    
    #CODE FOR ESN ERROR
    Wout <- train_Wout(y = train_data[, 2:number_of_train[i], drop = FALSE],
                       u = train_data[, 1:(number_of_train[i] - 1), drop = FALSE],
                       W = W, Win = Win, x = train_x[, 1:(number_of_train[i] - 1)],
                       alpha = alpha, beta = beta, print_det = print_det)
    
    esn_pred_val <- make_esn_predictions(u = test_data, W = W, Win = Win, x = test_x, alpha = alpha, Wout = Wout)
    MSE_array_esn[i, rho_idx] <- MSE_error(test_data[, 3:N_test, drop = FALSE], esn_pred_val[, 2:(N_test - 1), drop = FALSE])
  }
}
print(proc.time() - time0)