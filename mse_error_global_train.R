library(CVXR)
library(MASS)
wd <- getwd()
source(paste(wd, 'var_utils.R', sep = '/'))
source(paste(wd, 'esn_utils.R', sep = '/'))
#SIM PARAMETERS
N_sim <- 1
#DATA PARAMETERS
N <- 10
Ny <- 10
Nu <- 10
PHI_ro = 0.8
N_test <- 10000
sigma <- 0.1
#ESN PARAMETERS
Nx <- 1000
w_dens <- 1
w_ro <- 0.5
alpha <- 0.5
beta <- 0.01


number_of_train <- c(25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000, 750000)
MSE_array_var <- matrix(0, N_sim, length(number_of_train))
MSE_array_esn <- matrix(0, N_sim, length(number_of_train))
phi <- matrix(rnorm(1:N), N, 1)
PHI <- generate_PHI_normal(N = N, mu = 0, sigma = 1, ro = PHI_ro)
global_train_data <-generate_var1(number_of_train[length(number_of_train)], phi, PHI, diag(sigma, 10, 10)) #4.5 secs
test_data <- generate_var1(N_test, phi, PHI, diag(sigma, 10, 10))
time0 <- proc.time()
W <- as.matrix(generate_W(Nx = Nx, density = w_dens, sigma = 1, ro = w_ro))
Win <- generate_Win(Nx = Nx, Nu = Nu, mu = 0, sigma = 1)
train_x <- calculate_x_c(u = global_train_data[, 1:number_of_train[length(number_of_train)] - 1], W = W, Win = Win, alpha = alpha)
test_x <- calculate_x_c(u = test_data, W = W, Win = Win, alpha = alpha)

for (sim in 1:N_sim){
  print(sim)
  for (i in 1:length(number_of_train)){
    #CODE FOR VAR ERROR
    train_data <- global_train_data[, 1:number_of_train[i]]
    PHI_est <- estimate_parameters(train_data)
    var_pred_val <- make_var_predictions(test_data, PHI_est)
    MSE_array_var[sim, i] <- MSE_error(test_data[, 2:N_test], var_pred_val)
    #CODE FOR ESN ERROR
    Wout <- train_Wout(y = train_data[, 2:number_of_train[i]],
                       u = train_data[, 1:number_of_train[i] - 1],
                       W = W, Win = Win, x = train_x[, 1:number_of_train[i] - 1],
                       alpha = alpha, beta = beta)

    esn_pred_val <- make_esn_predictions(u = test_data, W = W, Win = Win, x = test_x, alpha = alpha, Wout = Wout)
    MSE_array_esn[sim, i] <- MSE_error(test_data[, 2:N_test], esn_pred_val[, 1:N_test - 1])
  }
}
print(proc.time() - time0)