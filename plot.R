library(rbokeh)
figure(width = 700, title = "Testing", 
       xlab="T", ylab="log-price", legend_location = "top_right") %>%
  ly_lines(number_of_train, MSE_array_esn[1, ], color = "red", width = 2, legend = "ESN") %>%
  ly_lines(number_of_train, MSE_array_var[1, ], color = "blue", width = 2, legend = "VAR") %>%
  x_axis(log = TRUE, label = 'SIZE OF TRAINING', )


#PARA PLOT
library(ggplot2)
library(reshape)
library(RColorBrewer)
it_array = beta_array
df <- as.data.frame(cbind(number_of_train, cbind(MSE_array_var, MSE_array_esn)))
df <- melt(df, id.vars = 'number_of_train')
gg <- ggplot(data=df,
       aes(x=number_of_train, y=value, colour=variable)) +
       theme_bw() + theme(text = element_text(size=14)) + theme(axis.text=element_text(size=rel(1.2))) +
       geom_line() + geom_point() + scale_x_log10() + xlab('SIZE OF TRAINING') + ylab('MSE ERROR') +
       scale_color_manual(name = 'BETA', labels = c('VAR model', it_array), values = brewer.pal(length(it_array) + 1, "Set1"))
plot(gg)
