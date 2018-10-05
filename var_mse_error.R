library(CVXR)
library(MASS)
wd <- getwd()
source(paste(wd, 'var_utils.R', sep = '/'))
source(paste(wd, 'esn_utils.R', sep = '/'))
N <- 10
Ny <- 10
Nu <- 10
N_test <- 10000
Nx <- 100
w_dens <- 0.2
w_ro <- 0.5
alpha <- 0.5
sigma <- 0.1

number_of_train <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000)
#number_of_train <- c(50000)
MSE_array <- matrix(0, 2, length(number_of_train))
time0 <- proc.time()
for (i in 1:length(number_of_train)){
  #CODE FOR VAR ERROR
  print(i)
  phi <- matrix(rnorm(1:N), N, 1)
  PHI <- generate_PHI_normal(N)
  train_data <-generate_var1(number_of_train[i], phi, PHI, diag(sigma, 10, 10))
  test_data <- generate_var1(N_test, phi, PHI, diag(sigma, 10, 10))
  PHI_est <- estimate_parameters(train_data)
  var_pred_val <- make_var_predictions(test_data, PHI_est)
  MSE_array[1, i] <- MSE_error(test_data[, 2:N_test], var_pred_val)
  #CODE FOR ESN ERROR
  W <- generate_W(Nx = Nx, density = w_dens, sigma = 1, ro = w_ro)
  Win <- generate_Win(Nx = Nx, Nu = Nu, mu = 0, sigma = 1)
  Wout <- train_Wout(y = train_data[, 2:number_of_train[i]],
                     u = train_data[, 1:number_of_train[i] - 1],
                     W = W, Win = Win, alpha = alpha)
  esn_pred_val <- make_esn_predictions(u = test_data, W = W, Win = Win, alpha = alpha, Wout = Wout)
  MSE_array[2, i] <- MSE_error(test_data[, 2:N_test], esn_pred_val[, 1:N_test - 1])
}
print(proc.time() - time0)

