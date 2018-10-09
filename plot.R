library(rbokeh)
figure(width = 700, title = "Testing", 
       xlab="T", ylab="log-price", legend_location = "top_right") %>%
  ly_lines(number_of_train, MSE_array_esn[1, ], color = "red", width = 2, legend = "ESN") %>%
  ly_lines(number_of_train, MSE_array_var[1, ], color = "blue", width = 2, legend = "VAR") %>%
  x_axis(log = TRUE, label = 'SIZE OF TRAINING')