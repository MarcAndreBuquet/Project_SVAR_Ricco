plot_figure_new <- function(data, data_lower = NULL, data_upper = NULL, horizon, nb_var, shock , var_names = var_names )
                
  
{
  IRFs_representation = as.list(numeric(length = nb_var))
  data_combined = as.list(numeric(length = nb_var))
  dim(IRFs_representation) = c(nb_var, 1)
  
  
  for (var in 1:nb_var) {data_combined[[var]] = cbind(data_lower[[var]], data_upper[[var]])}
      
  
  for (var in 1:nb_var) {
    
    temp_data = ggplot() + 
      geom_spline(data = as.data.frame(data_lower[[var]]) %>% setNames(., "IRFs")  , aes(x = 1:horizon, y = IRFs), linetype = "dotted", spar = 0.675) + 
      geom_spline(data = as.data.frame(data_upper[[var]]) %>% setNames(., "IRFs")  , aes(x = 1:horizon, y = IRFs), linetype = "dotted", spar = 0.675) 
    
    temp_data_1 <- ggplot_build(temp_data)
    temp_data_2 <- data.frame(x = temp_data_1$data[[1]]$x,
                              ymin = temp_data_1$data[[1]]$y,
                              ymax = temp_data_1$data[[2]]$y) 
    IRFs_representation[[var]] = ggplot() + 
      theme_bw() +
      geom_spline(data = as.data.frame(data_lower[[var]]) %>% setNames(., "IRFs")  , aes(x = 1:horizon, y = IRFs), linetype = "dotted", color = "black", spar = Pspar) + 
      geom_spline(data = as.data.frame(data_upper[[var]]) %>% setNames(., "IRFs")  , aes(x = 1:horizon, y = IRFs), linetype = "dotted", color = "black", spar = Pspar) + 
      labs(y = as.character(shock) , x = as.character( var_names[var]   )) +
      theme(axis.text.x= element_text(size=SL-2),axis.text.y= element_text(size=SL-2))+
      theme(axis.title.x = element_text(size=SL+3),axis.title.y= element_text(size=SL+3))+
      #geom_ribbon(data = as.data.frame(data_combined[var, shock]) %>% set_names(., "Lower_bound","Upper_bound"), aes(x=1:horizon, ymax = predict(loess(Upper_bound ~ (1:horizon))), ymin= predict(loess(Lower_bound ~ (1:horizon))), fill="pink", alpha=.2 )) +:
      geom_ribbon(data = temp_data_2, aes(x = x, ymin = ymin, ymax = ymax), fill = "pink", alpha = 0.2) + 
      geom_spline(data = as.data.frame(data[[var]]) %>% setNames(., "IRFs")  , aes(x = 1:horizon, y = IRFs), color = "black", spar = Pspar) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "red")
    
  }
  
  IRFs_final = gridExtra::grid.arrange(grobs = IRFs_representation, nrow = 1)

  
  
  return(IRFs_final)
  
  
}
