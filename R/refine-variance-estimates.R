refine_S_g_estimates <- function(S_g_list,
                                 Y_g_list,
                                 N_g_list,
                                 g_list,
                                 t_list){

  #This function takes a list of preliminary S_g estimates using the sample covariance for cohort g
  # It then updates the covariance estimates by pooling data across cohorts for all t<g

  # Create an estimate of the variance of untreated Y for all cohorts <= t
  create_pooled_S_for_t <- function(t){
    g_after_t_index <- which(g_list > t)
    time_leq_t_index <- which(t_list <= t)

    N_pooled <- purrr::reduce(N_g_list[g_after_t_index],
                              sum)
    #Ybar_pooled <- 1/N_pooled * reduce( map(.x = g_after_t_index, ~ N_g_list[[.x]] * Y_g_list[[.x]]) , sum)
    Ybar_pooled <- 1/N_pooled * Reduce(x = purrr::map2(.x = N_g_list,
                                                       .y = Y_g_list,
                                                       .f = ~ .x * .y[time_leq_t_index]),
                                       f= '+' )
    # Compute \sum_{i: g_i = g} (Y_i - Y_g) (Y_i - Y_g)'
    sum_squared_residuals_g_functions <- function(gIndex){
      S_g <- S_g_list[[gIndex]]
      Y_g <- Y_g_list[[gIndex]]
      N_g <- N_g_list[[gIndex]]

      S_g_sub <- S_g[time_leq_t_index, time_leq_t_index]
      Y_g_sub <- Y_g[time_leq_t_index]

      #This is sum of squared residuals using the bias-variance decomposition
      sum_squared_residuals_g <- (N_g - 1) * S_g_sub + N_g * (Y_g_sub - Ybar_pooled) %*% t(Y_g_sub - Ybar_pooled)

      return(sum_squared_residuals_g)
    }

    #Sum the residuals for all g > t, then divide by N_pooled
    pooled_S <- 1/(N_pooled - 1) * base::Reduce(x = purrr::map(.x = g_after_t_index,
                                                               .f = sum_squared_residuals_g_functions),
                                          sum)

    return(pooled_S)
  }

  refine_S_g <- function(gIndex){
    g <- g_list[gIndex]
    S_g <- S_g_list[[gIndex]]
    S_g_modified <- S_g

    #If no pre-treatment periods are in data, return S_g
    if(g <= min(t_list)){
      return(S_g)
      }

    #Create indiex of periods before g
    t_less_than_g_index <- which(t_list < g)

    #Start with S_g
    #For each period t < g, starting with t=g-1, replace S_g[ 1:t, 1:t ] with the pooled estimate for all cohorts treated after t
    for(tIndex in sort(t_less_than_g_index, decreasing = T) ){
      S_t_pooled <- create_pooled_S_for_t(t_list[tIndex])
      S_g_modified[1:tIndex, 1:tIndex] <- S_t_pooled}

    return(S_g_modified)
  }

 S_g_modified_list <- purrr::map(.x = 1:length(g_list),
                                 .f = refine_S_g)

 return(S_g_modified_list)
}


## tests ----


#refine_S_g_estimates(S_g_list, Ybar_g_list, N_g_list, g_list, t_list)
