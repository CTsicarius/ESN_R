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
N <- 10
Ny <- 10
Nu <- 10
PHI_ro = 0.2
N_test <- 10000
sigma <- 0.1
#ESN PARAMETERS
Nx <- 100
w_ro <- 0.001
w_dens <- 0.2
alpha <- 1
beta <- 150

set.seed(42)
number_of_train <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000)
#rho_array <- c(0, 0.0001, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1)
rho_array <- c(0, 0.01, 0.1, 0.25, 0.5, 0.75, 0.9, 1)

MSE_array_var <- matrix(0, length(number_of_train), 1)
MSE_array_esn <- matrix(0, length(number_of_train), length(rho_array))
phi <- matrix(rnorm(1:N), N, 1)
PHI1 <- generate_PHI_normal(N = N, mu = 0, sigma = 1, ro = PHI_ro)
PHI2 <- generate_PHI_normal(N = N, mu = 0, sigma = 1, ro = PHI_ro)
global_train_data <-generate_var2(number_of_train[length(number_of_train)], phi, PHI1, PHI2, diag(sigma, 10, 10)) #4.5 secs
test_data <- generate_var2(N_test, phi, PHI1, PHI2, diag(sigma, 10, 10))
time0 <- proc.time()
print('OK1')
Win <- generate_Win(Nx = Nx, Nu = Nu, mu = 0, sigma = 1)
for (rho_idx in 1:length(rho_array)){
  W <- as.matrix(generate_W(Nx = Nx, density = w_dens, sigma = 1, ro = rho_array[rho_idx]))
  train_x <- calculate_x_c(u = global_train_data[, 1:number_of_train[length(number_of_train)] - 1], W = W, Win = Win, alpha = alpha)
  test_x <- calculate_x_c(u = test_data, W = W, Win = Win, alpha = alpha)
  for (i in 1:length(number_of_train)){
    #CODE FOR VAR ERROR
    train_data <- global_train_data[, 1:number_of_train[i]]
    PHI_est <- estimate_var2_parameters(train_data)
    var_pred_val <- make_var2_predictions(test_data, PHI_est)
    MSE_array_var[i,  1] <- MSE_error(test_data[, 3:N_test], var_pred_val)
    #CODE FOR ESN ERROR
    Wout <- train_Wout(y = train_data[, 3:number_of_train[i]],
                       u = train_data[, 1:(number_of_train[i] - 2)],
                       W = W, Win = Win, x = train_x[, 1:(number_of_train[i] - 2)],
                       alpha = alpha, beta = beta_array[beta_idx], print_det = print_det)
    
    esn_pred_val <- make_esn_predictions(u = test_data, W = W, Win = Win, x = test_x, alpha = alpha, Wout = Wout)
    MSE_array_esn[i, rho_idx] <- MSE_error(test_data[, 3:N_test], esn_pred_val[, 1:(N_test - 2)])
  }
}
print(proc.time() - time0)