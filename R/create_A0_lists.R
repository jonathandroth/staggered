

create_A0_list_for_ATE_tg <- function(t,
                                      g,
                                      g_list,
                                      t_list,
                                      N_g_list,
                                      use_last_treated_only = FALSE){
  numPeriods <- length(t_list)
  #if(t < g){warning("t is less than g. ATE(t,g) is zero by assumption")}
  if(t >= max(g_list)){
    stop("t is greater than max(g)-1; ATE(t,g) is not identified.")
  }
  #Create A0s for ATT_{t,g}
  treated_cohort_index <- which(g_list == g)

  N_treated <- N_g_list[[treated_cohort_index]]

  #Create A0 for the treated units
  A0_g_treated_fn <- function(gIndex){
    A0_g <- matrix(0, nrow = 1, ncol = numPeriods)

    if(gIndex == treated_cohort_index){
      g <- g_list[gIndex]
      N_g <- N_g_list[[gIndex]]
      event_time_index <- which(t_list == g-1)
      A0_g[event_time_index] <- N_g / N_treated
    }
    return(A0_g)
  }


  if(!use_last_treated_only){
    #Create a list of which cohorts are eligible to be controls for each of the cohorts
    #This will be a null list if not eligible
    control_cohort_indices <- purrr::map(.x = 1:length(g_list),
                                         .f = ~ which(g_list > t ))
  }else{
    #If use_last_treated_only, compare only to the last treated cohort (i.e. max(G))
    control_cohort_indices <- purrr::map(.x = 1:length(g_list),
                                         .f = ~ which( (g_list > t) & (g_list ==  max(g_list) )) )
  }



  N_control_cohorts <- purrr::map_dbl(.x = control_cohort_indices,
                                      .f = ~ sum(unlist(N_g_list[.x])) )

  createControlWeights_helper <- function(control_g_index,
                                          treated_g_index = treated_cohort_index){
    #This function creates the weights for the cohort in control_g_index as a controls for the cohort in treated_g_index
    control_weights <- matrix(0,nrow =1, ncol = numPeriods)

    #If control_g is not a valid control, return 0
    if(! (control_g_index %in% control_cohort_indices[[treated_g_index]]) ){
      return(control_weights)
    }

    g_treated <- g_list[treated_g_index]
    N_g_control <- N_g_list[[control_g_index]]
    N_g_treated <- N_g_list[[treated_g_index]]

    t_control_index <- which(t_list == g-1)

    control_weights[t_control_index] <- -N_g_control / N_control_cohorts[treated_g_index]
    return(control_weights)
  }

  # createControlWeights <- function(control_g_index){
  #   #Create the control weights for control_g_index using all treated cohorts, then sum them
  #   controlWeights_g <- Reduce( x=  purrr::map(.x = eligible_cohort_index,
  #                                       ~createControlWeights_helper(control_g_index, .x)  ) ,
  #                               f= '+' )
  #   return(controlWeights_g)
  # }

  A0_treated_list <- purrr::map(.x = 1:length(g_list),
                                A0_g_treated_fn)
  A0_control_list <- purrr::map(.x = 1:length(g_list),
                                createControlWeights_helper )
  A0_list <- purrr::map2(.x = A0_treated_list,
                         .y = A0_control_list,
                         .f = ~ .x + .y)
  return(A0_list)
}


create_A0_list_for_event_study <- function(eventTime,
                                           g_list,
                                           t_list,
                                           N_g_list,
                                           use_last_treated_only = FALSE){

  #Create A0s for an ``event-study'' coefficient at lag eventTime
  # This estimand is the average treatment effects for units eventTime periods after first being treated
  # The average is taken over all cohorts g such that there is some untreated cohort at g+eventTime
  # Cohorts are weighted by the cohort size (N_g)
  # This function returns the pre-period differences used as control in CS
  # i.e, the estimate of the the difference in period g-1 btwn cohort g and units not-yet-treated at g+eventTime
  # If use_last_treated_only = TRUE, then only the last cohort is used as a control

  maxG <- max(g_list)
  eligible_cohort_index <- which( ((g_list + eventTime) < maxG ) & ((g_list + eventTime) <= max(t_list) ) )

  if(length(eligible_cohort_index) == 0){
    stop("There are no comparison cohorts for the given eventTime")
  }

  N_eligible <- Reduce(x = N_g_list[eligible_cohort_index], sum)

  A0_lists <- purrr::map(.x = eligible_cohort_index,
                         .f = ~ scalar_product_lists(N_g_list[[.x]]/N_eligible ,
                                                     create_A0_list_for_ATE_tg(t =g_list[.x]+eventTime,
                                                                               g = g_list[.x],
                                                                               g_list = g_list ,
                                                                               t_list = t_list,
                                                                               N_g_list = N_g_list,
                                                                               use_last_treated_only = use_last_treated_only) ) )

  if(length(eligible_cohort_index) == 1){
    A0_list <- A0_lists[[1]]
  }else{
    A0_list <- sum_of_lists(A0_lists)
  }
  return(A0_list)
}




create_A0_list_for_ATE_calendar_t <- function(t,
                                              g_list,
                                              t_list,
                                              N_g_list,
                                              use_last_treated_only = FALSE){

  treated_by_t_indices <- which(g_list <= t)
  N_total_treated <- sum( unlist(N_g_list[treated_by_t_indices]) )

  A0_lists <- purrr::map(.x = treated_by_t_indices,
                         .f = ~ scalar_product_lists(N_g_list[[.x]]/N_total_treated ,
                                                     create_A0_list_for_ATE_tg(t =t,
                                                                               g = g_list[.x],
                                                                               g_list = g_list ,
                                                                               t_list = t_list,
                                                                               N_g_list = N_g_list,
                                                                               use_last_treated_only = use_last_treated_only) ) )

  if(length(treated_by_t_indices) == 1){
    A0_list <- A0_lists[[1]]
  }else{
    A0_list <- sum_of_lists(A0_lists)
  }
  return(A0_list)
}




create_A0_list_for_ATE_cohort_g <- function(g,
                                            g_list,
                                            t_list,
                                            N_g_list,
                                            use_last_treated_only = FALSE){

  treated_period_indices <- which((t_list >= g) & (t_list < max(g_list)))
  T_treated <- length( t_list[treated_period_indices] )

  A0_lists <- purrr::map(.x = treated_period_indices,
                         .f = ~ scalar_product_lists(1/T_treated ,
                                                     create_A0_list_for_ATE_tg(t =t_list[.x],
                                                                               g = g,
                                                                               g_list = g_list ,
                                                                               t_list = t_list,
                                                                               N_g_list = N_g_list,
                                                                               use_last_treated_only = use_last_treated_only) ) )

  if(T_treated == 1){
    A0_list <- A0_lists[[1]]
  }else{
    A0_list <- sum_of_lists(A0_lists)
  }

  return(A0_list)
}




create_A0_list_for_cohort_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_list,
                                                  use_last_treated_only = FALSE){

  g_eligible_index <-  which(g_list < max(g_list) & g_list <= max(t_list))

  N_total_eligible <- sum(unlist(N_g_list[g_eligible_index]))

  A0_lists <- purrr::map(.x = g_eligible_index,
                         .f = ~ scalar_product_lists(as.numeric(N_g_list[[.x]])/N_total_eligible,
                                                     create_A0_list_for_ATE_cohort_g(g = g_list[.x],
                                                                                     g_list = g_list,
                                                                                     t_list = t_list,
                                                                                     N_g_list = N_g_list,
                                                                                     use_last_treated_only = use_last_treated_only) ) )


  A0_list <- sum_of_lists(A0_lists)

  return(A0_list)
}




create_A0_list_for_calendar_average_ATE <- function(g_list,
                                                    t_list,
                                                    N_g_list,
                                                    use_last_treated_only = FALSE){

  t_eligible_index <-  which(t_list >= min(g_list) & t_list < max(g_list))

  T_eligible <- length(t_eligible_index)

  A0_lists <- purrr::map(.x = t_eligible_index,
                         .f = ~ scalar_product_lists(1/T_eligible,
                                                     create_A0_list_for_ATE_calendar_t(t = t_list[.x],
                                                                                       g_list = g_list,
                                                                                       t_list = t_list,
                                                                                       N_g_list = N_g_list,
                                                                                       use_last_treated_only = use_last_treated_only) ) )


  A0_list <- sum_of_lists(A0_lists)

  return(A0_list)
}



create_A0_list_for_simple_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_list,
                                                  use_last_treated_only = FALSE){

  # AVOID NOTE ON CRAN
  g <- NULL

  #Create a df with all the (g,t) pairs for which ATE is identified
  gt_df <- purrr::cross_df( list(g = g_list,
                                 t = t_list) )
  gt_df <- gt_df %>% dplyr::filter(t >= g,
                                   t< max(g_list))

  #Join in N_g for each of these pairs
  gt_df <- dplyr::left_join( gt_df,
                      data.frame(g= g_list,
                                 N_g = unlist(N_g_list)),
                      by = "g")

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- sum(gt_df$N_g)

  A0_lists <- purrr::map(.x = 1:NROW(gt_df),
                         .f = ~ scalar_product_lists(gt_df$N_g[.x]/N_total,
                                                     create_A0_list_for_ATE_tg(t = gt_df$t[.x],
                                                                               g= gt_df$g[.x],
                                                                               g_list = g_list ,
                                                                               t_list = t_list,
                                                                               N_g_list = N_g_list,
                                                                               use_last_treated_only = use_last_treated_only) ) )


  A0_list <- sum_of_lists(A0_lists)

  return(A0_list)
}

