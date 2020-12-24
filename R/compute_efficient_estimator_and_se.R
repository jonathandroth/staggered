
#' @title Calculate group level summary statistics
#' @description This function computes the mean-vector and covariance matrix of the outcomes for each cohort, where a cohort g is a group of units first treated in period g
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @return Y_bar_list A list of the means of the outcomes for each cohort g
#' @return S_g_list A list of covariance matrices for the outcomes for each cohort g
#' @return N_g_list A list of the number of observations for each cohort g
#' @return g_list A list of when the cohorts were first treated
#' @return t_list A list of the the time periods for the outcome. The vector of outcomes corresponds with this order.
#' @export
compute_g_level_summaries <- function(df, refine_S_g = T){

  g_list <- sort(unique(df$g))
  t_list <- sort(unique(df$t))
  compute_Ybar_Sbar_g <- function(g){
    #Filter to observations in cohort g
    dfg <- df %>% dplyr::filter(g == !!g)

    #Reshape so that Y_{it} is a column for all t
    dfg <- dfg  %>%
      reshape2::dcast(i ~ t, value.var = "y") %>%
      dplyr::select(-i)

    #Order the columns ascending in time
    tVec <- as.numeric(colnames(dfg))
    dfg <- dfg[,order(tVec)]

    #Calculate the number of observation N_g
    N_g <- NROW(dfg)

    #Convert to matrix for taking a covariance
    dfg <- dfg %>% as.matrix()

    #Compute means Ybar_g and covariance S_g
    Ybar_g <- colMeans(dfg)
    S_g <- var(dfg)


    return(list(Ybar_g = Ybar_g, S_g = S_g, N_g = N_g ) )
  }

  resultsList <- purrr::map(.x = g_list, .f = compute_Ybar_Sbar_g)

  Ybar_g_List <- purrr::map(.x = resultsList, .f = ~.x$Ybar_g)
  S_g_List <- purrr::map(.x = resultsList, .f = ~.x$S_g)
  N_g_List <- purrr::map(.x = resultsList, .f = ~.x$N_g)

  if(refine_S_g){S_g_List <- refine_S_g_estimates(S_g_list = S_g_List, Y_g_list = Ybar_g_List, N_g_list = N_g_List, g_list = g_list, t_list = t_list)}

  return( list( Ybar_g_List = Ybar_g_List, S_g_List = S_g_List, N_g_List = N_g_List, g_list = g_list, t_list = t_list ) )
}


compute_Thetahat0 <- function(Ybar_g_list, A_theta_list){
  A_theta_Ybar_list <- purrr::map2(.x = Ybar_g_list, .y = A_theta_list, .f = ~.y %*% .x)
  Thetahat0 <- purrr::reduce(.x = A_theta_Ybar_list, .f = sum)
  return(Thetahat0)
}

compute_Xhat <- function(Ybar_g_list, A_0_list){
  A_0_Ybar_list <- purrr::map2(.x = Ybar_g_list, .y = A_0_list, .f = ~.y %*% .x)
  Xhat <- purrr::reduce(.x = A_0_Ybar_list, .f = sum)
  return(Xhat)
}


compute_Betastar <- function(Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list = NULL){

  if(is.null(Xvar_list)){
    #Xvar_list <- pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*% S %*% t(A0) ) } )
    Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , t(A0) ) ) } )
  }

  Xvar <- purrr::reduce(.x = Xvar_list, .f = sum)

  X_theta_cov_list <- purrr::pmap(.l = list(A_0_list, A_theta_list, S_g_list, N_g_list) , .f = function(A0,A_theta,S,N){ return(1/N * A0 %*% S %*% t(A_theta) ) } )
  X_theta_cov <- purrr::reduce(.x = X_theta_cov_list, .f = sum)

  betastar <- solve(Xvar) %*% X_theta_cov
  return(betastar)
}


compute_Thetahat_beta <- function(beta, Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list = NULL){
  thetahat0 <- compute_Thetahat0(Ybar_g_list, A_theta_list)
  Xhat <- compute_Xhat(Ybar_g_list, A_0_list)
  thetahatbeta <- thetahat0 - t(Xhat) * beta

  return(thetahatbeta)
}

compute_Thetahat_Star <- function(Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list){

  betastar <- compute_Betastar(Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list)

  thetahatstar <- compute_Thetahat_beta(beta = betastar, Ybar_g_list = Ybar_g_list, A_theta_list = A_theta_list, A_0_list = A_0_list, S_g_list = S_g_list, N_g_list = N_g_list)

  return(thetahatstar)
}

compute_se_Thetahat_beta_conservative <- function(beta, Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list,Xvar_list = NULL, ...){

  if(is.null(Xvar_list)){
    #Xvar_list <- pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*% S %*% t(A0) ) } )
    Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , t(A0) ) ) } )
  }
  Xvar <- purrr::reduce(.x = Xvar_list, .f = sum)

  X_theta_cov_list <- purrr::pmap(.l = list(A_0_list, A_theta_list, S_g_list, N_g_list) , .f = function(A0,A_theta,S,N){ return(1/N * A0 %*% S %*% t(A_theta) ) } )
  X_theta_cov <- purrr::reduce(.x = X_theta_cov_list, .f = sum)

  thetaVar_conservative_list <- purrr::pmap(.l = list(A_theta_list, S_g_list, N_g_list) , .f = function(A_theta,S,N){ return(1/N * A_theta %*% S %*% t(A_theta) ) } )
  thetaVar_conservative <- purrr::reduce(.x = thetaVar_conservative_list, .f = sum)

  varhat_conservative <- thetaVar_conservative + t(beta) %*% Xvar %*% beta - 2* X_theta_cov %*% beta
  se_conservative <- sqrt(varhat_conservative)
  return(se_conservative)
}

computeGMin <- function(A_theta_list,g_list){
  A_theta_is_nonzero <- purrr::map_lgl(.x = A_theta_list, .f = ~max( abs(.x)) != 0 )
  min_nonzero_index <- min( which(A_theta_is_nonzero) )
  g_min <- g_list[min_nonzero_index]

  if(length(g_min) == 0){ g_min = g_list[1] }

  return(g_min)
}



compute_se_Thetahat_beta <- function(beta, Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, g_list, t_list, Xvar_list = NULL){

  seConservative <- compute_se_Thetahat_beta_conservative(beta, Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list =  Xvar_list)

  gMin <- computeGMin(A_theta_list = A_theta_list, g_list = g_list)
  g_geq_gMin_index <- which(g_list >= gMin) #which indices of g are >= gMin

  tMin <- min(t_list)
  tMax <- max(t_list)


  #Right now, we're assuming first period is 1. May want to relax this
  if(gMin == tMin){return(seConservative)}

  #Create matrix M that selects the rows of S_g correspondign with t< g_min
  M <- matrix(0, nrow = gMin -tMin, ncol = NCOL(S_g_list[[1]]) )
  diag(M) <- 1



  betahat_g_list <- purrr::pmap( .l = list( A_theta_list[g_geq_gMin_index], A_0_list[g_geq_gMin_index], S_g_list[g_geq_gMin_index]  ),
                          .f = function(A_theta,A_0,S_g){ MASS::ginv(M %*% S_g %*% t(M)) %*% M %*% S_g %*% t(A_theta) } )

  betahat_g_sum <- purrr::reduce(betahat_g_list, sum)

  #Compute the average of M S_g M' for all g <= t_min
  avg_MSM_list <- purrr::map(.x = S_g_list[g_geq_gMin_index], .f = function(S){ M %*% S %*% t(M) } )
  avg_MSM <- purrr::reduce(avg_MSM_list, sum) / length(avg_MSM_list)

  #calculate the adjustment factor for the conservative variance
  N <- purrr::reduce(N_g_list, sum)
  adjustmentFactor <- 1/N * t(betahat_g_sum) %*% avg_MSM %*% betahat_g_sum



  var_conservative <- seConservative^2
  if(var_conservative - adjustmentFactor < 0){
    warning("var_conservative is less than adjustmentFactor")
    se_adjusted <- 0
  }else{
    se_adjusted <- sqrt( var_conservative - adjustmentFactor )
  }
  return(se_adjusted)
}


create_A0_list <- function(g_list, t_list){

  createAtilde0_g <- function(g_index){
    g <- g_list[g_index]
    t_less_than_g_index <- which(t_list < g)

    #If no periods less than g, there are no restrictions, return null
    if(length(t_less_than_g_index) == 0){
      return(NULL)
    }

    #Otherwise, Atilde is a matrix with rows= # periods before g and T columns.
    # It has 1s on the diagonal and 0s elsewhere. I.e. Atilde selects rows of Y corresponding with t<g
    Atilde0_g <- matrix(0, nrow = length(t_less_than_g_index), ncol = length(t_list))
    diag(Atilde0_g) <- 1

    return(Atilde0_g)
  }

  Atilde0_list <- purrr::map(.x = 1:(length(g_list)-1), .f = createAtilde0_g)

  A0_list <- purrr::map(.x = 1:length(Atilde0_list),
                 .f = function(i){
                   #A0_g is a block matrix that stacks the Atilde_g's (for g < gmax) and multiplies by 0 for all g' \neq g
                   Atilde_times_indicator_list <- purrr::map(.x = 1:length(Atilde0_list), .f = ~(i == .x) * Atilde0_list[[.x]]  )
                   A0 <- purrr::reduce(Atilde_times_indicator_list, rbind)
                   return(A0)
                 }  )
  if(length(A0_list) == 1){
    A0_gmax <- -A0_list[[1]]}
  else{
    A0_gmax <- -Reduce(x = A0_list, f = '+')
  }

  A0_list[[length(A0_list) + 1]] <- A0_gmax
  return(A0_list)
}


create_Atheta_list_for_event_study <- function(eventTime, g_list, t_list, N_g_list){

  #Create A_thetas for an ``event-study'' coefficient at lag eventTime
  # This is the average treatment effects for units eventTime periods after first being treated
  # The average is taken over all cohorts g such that there is some untreated cohort at g+eventTime
  # Cohorts are weighted by the cohort size (N_g)

  maxG <- max(g_list)
  eligible_cohort_index <- which(g_list + eventTime < maxG)

  if(length(eligible_cohort_index) == 0){stop("There are no comparison cohorts for the given eventTime")}

  N_eligible <- Reduce(x = N_g_list[eligible_cohort_index], sum)

  numPeriods <- length(t_list)

  #Create A_theta for the treated units
  A_theta_g_treated_fn <- function(gIndex){
    A_theta_g <- matrix(0, nrow = 1, ncol = numPeriods)

    if(gIndex %in% eligible_cohort_index){
      g <- g_list[gIndex]
      N_g <- N_g_list[[gIndex]]
      event_time_index <- which(t_list == g + eventTime)
      A_theta_g[event_time_index] <- N_g / N_eligible
    }
    return(A_theta_g)
  }


  #Create a list of which cohorts are eligible to be controls for each of the cohorts
  #This will be a null list if not eligible
  control_cohort_indices <- purrr::map(.x = 1:length(g_list),
                                .f = ~ which(g_list > g_list[[.x]] + eventTime ))

  N_control_cohorts <- purrr::map_dbl(.x = control_cohort_indices, .f = ~ sum(unlist(N_g_list[.x])) )

  createControlWeights_helper <- function(control_g_index, treated_g_index){
    #This function creates the weights for hte cohort in control_g_index as a controls for the cohort in treated_g_index

    control_weights <- matrix(0,nrow =1, ncol = numPeriods)

    #If control_g is not a valid control, return 0
    if(! (control_g_index %in% control_cohort_indices[[treated_g_index]]) ){return(control_weights)}

    g_treated <- g_list[treated_g_index]
    N_g_control <- N_g_list[[control_g_index]]
    N_g_treated <- N_g_list[[treated_g_index]]

    t_control_index <- which(t_list == g_treated + eventTime)

    control_weights[t_control_index] <- -N_g_control / N_control_cohorts[treated_g_index] * N_g_treated / N_eligible
    return(control_weights)
  }

  createControlWeights <- function(control_g_index){
    #Create the control weights for control_g_index using all treated cohorts, then sum them
    controlWeights_g <- Reduce( x=  purrr::map(.x = eligible_cohort_index,
                                        ~createControlWeights_helper(control_g_index, .x)  ) ,
                                f= '+' )
    return(controlWeights_g)
  }

  A_theta_treated_list <- purrr::map(.x = 1:length(g_list), A_theta_g_treated_fn)
  A_theta_control_list <- purrr::map(.x = 1:length(g_list), createControlWeights )
  A_theta_list <- purrr::map2(.x = A_theta_treated_list, .y = A_theta_control_list, .f = ~ .x + .y)
  return(A_theta_list)
}


sum_of_lists <- function(.l){
  #Given a list of lists all of the same length, this sums all the elements
  if(length(.l) == 1){return(.l)}
  results <- purrr::reduce(.x = .l, .f = function(l1,l2){purrr::map2(l1,l2, .f = ~.x + .y)} )
  return(results)
}

scalar_product_lists <- function(c,.l){
  #Takes a constant c and list .l, returns a list with c times each element of l
  results <- purrr::map(.x = .l, .f = ~c * .x)
  return(results)
}


create_Atheta_list_for_ATE_tg <- function(t, g, g_list, t_list, N_g_list){
  numPeriods <- length(t_list)
  if(t < g){warning("t is less than g. ATE(t,g) is zero by assumption")}
  if(t >= max(g_list)){stop("t is greater than max(g)-1; ATE(t,g) is not identified.")}
  #Create A_thetas for ATT_{t,g}
  treated_cohort_index <- which(g_list == g)

  N_treated <- N_g_list[[treated_cohort_index]]

  #Create A_theta for the treated units
  A_theta_g_treated_fn <- function(gIndex){
    A_theta_g <- matrix(0, nrow = 1, ncol = numPeriods)

    if(gIndex == treated_cohort_index){
      g <- g_list[gIndex]
      N_g <- N_g_list[[gIndex]]
      event_time_index <- which(t_list == t)
      A_theta_g[event_time_index] <- N_g / N_treated
    }
    return(A_theta_g)
  }


  #Create a list of which cohorts are eligible to be controls for each of the cohorts
  #This will be a null list if not eligible
  control_cohort_indices <- purrr::map(.x = 1:length(g_list),
                                .f = ~ which(g_list > t ))

  N_control_cohorts <- purrr::map_dbl(.x = control_cohort_indices, .f = ~ sum(unlist(N_g_list[.x])) )

  createControlWeights_helper <- function(control_g_index, treated_g_index = treated_cohort_index){
    #This function creates the weights for the cohort in control_g_index as a controls for the cohort in treated_g_index
    control_weights <- matrix(0,nrow =1, ncol = numPeriods)

    #If control_g is not a valid control, return 0
    if(! (control_g_index %in% control_cohort_indices[[treated_g_index]]) ){return(control_weights)}

    g_treated <- g_list[treated_g_index]
    N_g_control <- N_g_list[[control_g_index]]
    N_g_treated <- N_g_list[[treated_g_index]]

    t_control_index <- which(t_list == t)

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

  A_theta_treated_list <- purrr::map(.x = 1:length(g_list), A_theta_g_treated_fn)
  A_theta_control_list <- purrr::map(.x = 1:length(g_list), createControlWeights_helper )
  A_theta_list <- purrr::map2(.x = A_theta_treated_list, .y = A_theta_control_list, .f = ~ .x + .y)
  return(A_theta_list)
}


create_Atheta_list_for_ATE_calendar_t <- function(t, g_list, t_list, N_g_list){

  treated_by_t_indices <- which(g_list <= t)
  N_total_treated <- sum( unlist(N_g_list[treated_by_t_indices]) )

  A_theta_lists <- purrr::map(.x = treated_by_t_indices,
                       .f = ~ scalar_product_lists(N_g_list[[.x]]/N_total_treated ,
                                                   create_Atheta_list_for_ATE_tg(t =t, g = g_list[.x], g_list = g_list ,t_list = t_list, N_g_list = N_g_list) ) )

  if(length(treated_by_t_indices) == 1){
    A_theta_list <- A_theta_lists[[1]]
  }else{
    A_theta_list <- sum_of_lists(A_theta_lists)
  }
  return(A_theta_list)
}




create_Atheta_list_for_ATE_cohort_g <- function(g, g_list, t_list, N_g_list){

  treated_period_indices <- which(t_list >= g & t_list < max(g_list))
  T_treated <- length( t_list[treated_period_indices] )

  A_theta_lists <- purrr::map(.x = treated_period_indices,
                       .f = ~ scalar_product_lists(1/T_treated ,
                                                   create_Atheta_list_for_ATE_tg(t =t_list[.x], g = g, g_list = g_list ,t_list = t_list, N_g_list = N_g_list) ) )

  if(T_treated == 1){
    A_theta_list <- A_theta_lists[[1]]
  }else{
    A_theta_list <- sum_of_lists(A_theta_lists)
  }

  return(A_theta_list)
}




create_Atheta_list_for_cohort_average_ATE <- function(g_list, t_list, N_g_list){

  g_eligible_index <-  which(g_list < max(g_list) & g_list <= max(t_list))

  N_total_eligible <- sum(unlist(N_g_list[g_eligible_index]))

  A_theta_lists <- purrr::map(.x = g_eligible_index,
                       .f = ~ scalar_product_lists(as.numeric(N_g_list[[.x]])/N_total_eligible,
                                                   create_Atheta_list_for_ATE_cohort_g(g = g_list[.x], g_list = g_list ,t_list = t_list, N_g_list = N_g_list) ) )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}




create_Atheta_list_for_calendar_average_ATE <- function(g_list, t_list, N_g_list){

  t_eligible_index <-  which(t_list >= min(g_list) & t_list < max(g_list))

  T_eligible <- length(t_eligible_index)

  A_theta_lists <- purrr::map(.x = t_eligible_index,
                       .f = ~ scalar_product_lists(1/T_eligible,
                                                   create_Atheta_list_for_ATE_calendar_t(t = t_list[.x], g_list = g_list ,t_list = t_list, N_g_list = N_g_list) ) )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}



create_Atheta_list_for_simple_average_ATE <- function(g_list, t_list, N_g_list){

  #Create a df with all the (g,t) pairs for which ATE is identified
  gt_df <- purrr::cross_df( list(g = g_list, t = t_list) )
  gt_df <- gt_df %>% filter(t >= g, t< max(g_list))

  #Join in N_g for each of these pairs
  gt_df <- left_join( gt_df, data.frame(g= g_list, N_g = unlist(N_g_list)), by = "g")

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- sum(gt_df$N_g)

  A_theta_lists <- purrr::map(.x = 1:NROW(gt_df),
                       .f = ~ scalar_product_lists(gt_df$N_g[.x]/N_total,
                                                   create_Atheta_list_for_ATE_tg(t = gt_df$t[.x], g= gt_df$g[.x], g_list = g_list ,t_list = t_list, N_g_list = N_g_list) ) )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}

#source(here("Code/refine-variance-estimates.R"))



#' @export
#' @useDynLib staggered
#' @title Calculate the efficient adjusted estimator in staggered rollout designs
#' @description This functions calculates the efficient estimator for staggered rollout designs proposed by Roth and Sant'Anna.
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param estimand The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. The parameter can be left blank if a custom parameter is provided in A_theta_list
#' @param A_theta_list This parameter allows for specifying a custom estimand, and should be left as NULL if estimand is specified. It is a list of matrices A_theta_g so that the parameter of interest is \sum_g A_theta_g Ybar_g, where Ybar_g = 1/N \sum_i Y_i(g)
#' @param beta Optional. A coefficient to use for covariate adjustment. If not specified, the plug-in optimal coefficient is used
#' @param betaType. An optional parameter describing the type of covariate adjustment used. Defaults to "betaStar" if the efficient estimator is used.
#' @return df A data.frame containing: thetahat (the point estimate), se (the standard error), se_conservative (the Neyman standard error), and betaType
calculate_adjusted_estimator_and_se <- function(df, estimand =NULL, A_theta_list = NULL, A_0_list = NA, beta = NA, betaType = ifelse( is.na(beta), "betaStar", as.character(beta)), refine_S_g = T){

  g_level_summaries <- compute_g_level_summaries(df, refine_S_g = refine_S_g)
  Ybar_g_list <- g_level_summaries$Ybar_g_List
  S_g_list <- g_level_summaries$S_g_List
  N_g_list <- g_level_summaries$N_g_List
  g_list <- g_level_summaries$g_list
  t_list <- g_level_summaries$t_list

  #If estimand is provided, calculate the appropriate A_theta_list
  if(!is.null(estimand)){
    if(estimand == "simple"){
      A_theta_list <- create_Atheta_list_for_simple_average_ATE(g_list = g_list, t_list = t_list,N_g_list = N_g_list)
    }else if(estimand == "cohort"){
      A_theta_list <- create_Atheta_list_for_cohort_average_ATE(g_list = g_list, t_list = t_list,N_g_list = N_g_list)
    }else if(estimand == "calendar"){
      A_theta_list <- create_Atheta_list_for_calendar_average_ATE(g_list = g_list, t_list = t_list,N_g_list = N_g_list)
    }
  }
  #If no valid estimand is provided and no A_theta_list, throw and error
  if(is.null(A_theta_list)){stop("Estimand must be one of simple, cohort, or calendar; or custom A_theta_list must be provided")}

  #Create A_0_list if a custom A_0_list is not provided
  if(is.na(A_0_list)){
    A_0_list <- create_A0_list(g_list = g_list, t_list = t_list)
  }

  Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , t(A0) ) ) } )
#  Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*%S %*% t(A0) )  } )

  if(is.na(beta)){
    beta <- compute_Betastar(Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list = Xvar_list)
  }


  thetahat <- compute_Thetahat_beta(beta = beta,Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list = Xvar_list)
  se_conservative <- compute_se_Thetahat_beta_conservative(beta = beta,Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, Xvar_list = Xvar_list)
  se <- compute_se_Thetahat_beta(beta = beta,Ybar_g_list, A_theta_list, A_0_list, S_g_list, N_g_list, g_list, t_list, Xvar_list = Xvar_list)

  df <- data.frame(thetahat = thetahat, se = se, se_conservative = se_conservative, betaType = betaType)

  return(df)
}
