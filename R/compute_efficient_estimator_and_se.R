
#' @title Calculate group level summary statistics
#' @description This function computes the mean-vector and covariance matrix of the outcomes for each cohort, where a cohort g is a group of units first treated in period g
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param is_balanced If true, the df has previously been balanced so this does not need to be done internally.
#' @return Y_bar_list A list of the means of the outcomes for each cohort g
#' @return S_g_list A list of covariance matrices for the outcomes for each cohort g
#' @return N_g_list A list of the number of observations for each cohort g
#' @return g_list A list of when the cohorts were first treated
#' @return t_list A list of the the time periods for the outcome. The vector of outcomes corresponds with this order.

compute_g_level_summaries <- function(df, is_balanced = TRUE){

  # Avoid NOTES in CRAN
  i <- NULL

  #Balance the panel (and throw a warning if original panel is unbalanced)
  if(!is_balanced){
    df <- balance_df(df)
  }
  g_list <- sort(unique(df$g))
  t_list <- sort(unique(df$t))

  #Reshape so that Y_{it} is a column for all t
  df <- df  %>%
    tidyr::pivot_wider(id_cols = c("i","g"), names_from = "t", values_from = "y") %>%
    dplyr::ungroup() %>%
    dplyr::select(-i)

  compute_Ybar_Sbar_g <- function(g){
    #Filter to observations in cohort g
    dfg <- df %>%
           dplyr::filter(g == !!g) %>%
           dplyr::select(-g)

    #Order the columns ascending in time
    tVec <- as.numeric(colnames(dfg))
    dfg <- dfg[,order(tVec)]

    #Calculate the number of observation N_g
    N_g <- base::NROW(dfg)

    #Convert to matrix for taking a covariance
    dfg <- dfg %>% as.matrix()

    #Compute means Ybar_g and covariance S_g
    Ybar_g <- base::colMeans(dfg)
    S_g <- coop::covar(dfg)
    # If covar is NA or NaN, replace it to 0
    #S_g[is.na(S_g)] <- 0


    return(list(Ybar_g = Ybar_g,
                S_g = S_g,
                N_g = N_g))
  }


  resultsList <- purrr::map(.x = g_list,
                            .f = compute_Ybar_Sbar_g)

  Ybar_g_List <- purrr::map(.x = resultsList,
                            .f = ~.x$Ybar_g)

  S_g_List <- purrr::map(.x = resultsList,
                         .f = ~.x$S_g)

  N_g_List <- purrr::map(.x = resultsList,
                         .f = ~.x$N_g)


  return( list( Ybar_g_List = Ybar_g_List,
                S_g_List = S_g_List,
                N_g_List = N_g_List,
                g_list = g_list,
                t_list = t_list ))
}



balance_df <- function(df){

  ## This function creates a balanced panel as needed for our analysis

  # It first checks if rows of the data are uniquely characterized by (i,t)
  # If there are multiple observations per (i,t), it throws an error
  # It also removes observations with missing y

  #It then removes observations i for which data is not available for all t

  # Avoid NOTES in CRAN
  i <- NULL
  y <- NULL
  numPeriods_i <- NULL
  t <- NULL



  numPeriods <- length(unique(df$t))


  ##Check that (i,t) is a unique identifier
  if(anyDuplicated(df[c("i","t")]) > 0 ){
    stop("There are multiple observations with the same (i,t) values. The panel should have a unique outcome for each (i,t) value.")
  }

  #Check if there are missign values for y, and remove them if so
  if(max(is.na(df$y)) > 0 ){
    df <- df %>% dplyr::filter(!is.na(y))
  }

  #Check if panel is balanced. If not, drop the unbalanced observations and throw a warning
  df <- df %>%
    dplyr::group_by(i) %>%
    dplyr::mutate(numPeriods_i = length(unique(t)))

  if(min(df$numPeriods_i) == numPeriods){
    return(df)
  }else{
    warning("Panel is unbalanced (or has missing values for some observations). Dropping observations with missing values of Y_{it} for some time periods. If you wish to include these observations, provide staggered with a df with imputed outcomes.")
    df <- df %>% dplyr::filter(numPeriods_i == numPeriods) %>% dplyr::select(-numPeriods_i)
    return(df)
  }
}

compute_Thetahat0 <- function(Ybar_g_list,
                              A_theta_list){

  A_theta_Ybar_list <- purrr::map2(.x = Ybar_g_list,
                                   .y = A_theta_list,
                                   .f = ~.y %*% .x)
  Thetahat0 <- purrr::reduce(.x = A_theta_Ybar_list,
                             .f = sum)
  return(Thetahat0)
}


#' @title Compute Xhat of pre-treatment differences
#' @description \code{compute_Xhat} computes the vector Xhat of pre-treatment differences given the list of cohort means
#' Ybar_g_list and the list of matrices A_0_list
#' @param Ybar_g_list Ybar_g_list
#' @param A_0_list A_0_list
#' @return Xhat the vector Xhat of pre-treatment differences to be used as regressors
compute_Xhat <- function(Ybar_g_list, A_0_list){
  A_0_Ybar_list <- purrr::map2(.x = Ybar_g_list,
                               .y = A_0_list,
                               .f = ~.y %*% .x)
  Xhat <- base::Reduce(x = A_0_Ybar_list,
                       f = '+')
  return(Xhat)
}


#' @title  Plug-in efficient Beta hat
#' @description \code{compute_Betastar} computes the plug-in efficient betahat
#' @param Ybar_g_list Ybar_g_list
#' @param A_theta_list A_theta_list
#' @param A_0_list A_0_list
#' @param S_g_list S_g_list
#' @param N_g_list N_g_list
#' @param Xvar_list Xvar_list
#' @return betastar Vector of plug-in efficient betahat estimates.


compute_Betastar <- function(Ybar_g_list,
                             A_theta_list,
                             A_0_list,
                             S_g_list,
                             N_g_list,
                             Xvar_list = NULL){

  if(is.null(Xvar_list)){
    #Xvar_list <- pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*% S %*% base::t(A0) ) } )
    Xvar_list <- purrr::pmap(.l = list(A_0_list,
                                       S_g_list, N_g_list),
                             .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , base::t(A0) ) ) }
    )
  }

  #Xvar <- purrr::reduce(.x = Xvar_list, .f = sum)
  Xvar <- base::Reduce(f = '+',
                       x= Xvar_list)

  X_theta_cov_list <- purrr::pmap(.l = list(A_0_list,
                                            A_theta_list,
                                            S_g_list,
                                            N_g_list),
                                  .f = function(A0,A_theta,S,N){ return(1/N * A0 %*% S %*% base::t(A_theta) ) }
  )
  #X_theta_cov <- purrr::reduce(.x = X_theta_cov_list, .f = sum)
  X_theta_cov <- base::Reduce(x = X_theta_cov_list, f = '+')

  #betastar <- solve(Xvar) %*% X_theta_cov
  #betastar <- MASS::ginv(Xvar) %*% X_theta_cov
  #betastar <- solve_least_squares_svd(Xvar,X_theta_cov)
  betastar <- solve_least_squares_normal(Xvar,X_theta_cov) #fast method of solving (Xvar)^-1 X_theta_cov
  return(betastar)
}


compute_Thetahat_beta <- function(beta,
                                  Ybar_g_list,
                                  A_theta_list,
                                  A_0_list,
                                  S_g_list,
                                  N_g_list,
                                  Xvar_list = NULL){

  thetahat0 <- compute_Thetahat0(Ybar_g_list,
                                 A_theta_list)

  Xhat <- compute_Xhat(Ybar_g_list,
                       A_0_list)

  thetahatbeta <- thetahat0 - base::crossprod(Xhat, beta)

  return(thetahatbeta)
}

compute_Thetahat_Star <- function(Ybar_g_list,
                                  A_theta_list,
                                  A_0_list,
                                  S_g_list,
                                  N_g_list){

  betastar <- compute_Betastar(Ybar_g_list,
                               A_theta_list,
                               A_0_list,
                               S_g_list,
                               N_g_list)

  thetahatstar <- compute_Thetahat_beta(beta = betastar,
                                        Ybar_g_list = Ybar_g_list,
                                        A_theta_list = A_theta_list,
                                        A_0_list = A_0_list,
                                        S_g_list = S_g_list,
                                        N_g_list = N_g_list)

  return(thetahatstar)
}

compute_se_Thetahat_beta_conservative <- function(beta,
                                                  Ybar_g_list,
                                                  A_theta_list,
                                                  A_0_list,
                                                  S_g_list,
                                                  N_g_list,
                                                  Xvar_list = NULL,
                                                  ...){

  if(is.null(Xvar_list)){
    #Xvar_list <- pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*% S %*% base::t(A0) ) } )
    Xvar_list <- purrr::pmap(.l = list(A_0_list,
                                       S_g_list,
                                       N_g_list),
                             .f = function(A0,S,N){
                               return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , base::t(A0) ) )
                             }
    )
  }
  Xvar <- base::Reduce(x = Xvar_list,
                       f = '+')

  X_theta_cov_list <- purrr::pmap(.l = list(A_0_list,
                                            A_theta_list,
                                            S_g_list,
                                            N_g_list),
                                  .f = function(A0,A_theta,S,N){
                                    return(1/N * A0 %*% S %*% base::t(A_theta) )
                                  }
  )
  X_theta_cov <- base::Reduce(x = X_theta_cov_list,
                              f = '+')

  thetaVar_conservative_list <- purrr::pmap(.l = list(A_theta_list,
                                                      S_g_list,
                                                      N_g_list)
                                            , .f = function(A_theta,S,N){
                                              return(1/N * A_theta %*% S %*% base::t(A_theta) )
                                            }
  )
  thetaVar_conservative <- base::Reduce(x = thetaVar_conservative_list,
                                        f = '+')

  varhat_conservative <- thetaVar_conservative + base::t(beta) %*% Xvar %*% beta - 2 * base::t(X_theta_cov) %*% beta

  if(varhat_conservative <0){
    warning("Calculated variance is less than 0. Setting to 0.")
    se_neyman <- 0
  }else{
    se_neyman <- sqrt(varhat_conservative)
  }
  return(se_neyman)
}

computeGMin <- function(A_theta_list,
                        g_list){

  A_theta_is_nonzero <- purrr::map_lgl(.x = A_theta_list,
                                       .f = ~max( abs(.x)) != 0 )

  min_nonzero_index <- min( which(A_theta_is_nonzero) )
  g_min <- g_list[min_nonzero_index]

  if(length(g_min) == 0){
    g_min = g_list[1]
  }

  return(g_min)
}



compute_se_Thetahat_beta <- function(beta,
                                     Ybar_g_list,
                                     A_theta_list,
                                     A_0_list,
                                     S_g_list,
                                     N_g_list,
                                     g_list,
                                     t_list,
                                     Xvar_list = NULL,
                                     return_beta_sum = FALSE){

  #This function computes the standard error, using the version sigma_** that adjust for pre-treatment covariates
  # If return_beta_sum = TRUE, it returns a list wit

  seConservative <- compute_se_Thetahat_beta_conservative(beta,
                                                          Ybar_g_list,
                                                          A_theta_list,
                                                          A_0_list,
                                                          S_g_list,
                                                          N_g_list,
                                                          Xvar_list =  Xvar_list)

  gMin <- computeGMin(A_theta_list = A_theta_list,
                      g_list = g_list)

  g_geq_gMin_index <- which(g_list >= gMin) #which indices of g are >= gMin

  tMin <- min(t_list)
  tMax <- max(t_list)


  #Right now, we're assuming first period is 1. May want to relax this
  if(gMin == tMin){
    if(!return_beta_sum){
      return(seConservative)
    }else{
      resultsList <- list(se = seConservative,
                          betahat_g_sum = 0,
                          avg_MSM = 0,
                          N= purrr::reduce(N_g_list,
                                           sum))
      return(resultsList)
    }
  }

  #Create matrix M that selects the rows of S_g correspondign with t< g_min
  M <- matrix(0, nrow = gMin -tMin,
              ncol = NCOL(S_g_list[[1]]) )
  diag(M) <- 1



  betahat_g_list <- purrr::pmap( .l = list( A_theta_list[g_geq_gMin_index],
                                            A_0_list[g_geq_gMin_index],
                                            S_g_list[g_geq_gMin_index]  ),
                                 .f = function(A_theta,A_0,S_g){
                                   MASS::ginv(M %*% S_g %*% base::t(M)) %*% M %*% S_g %*% base::t(A_theta)
                                 }
  )

  betahat_g_sum <- base::Reduce(x = betahat_g_list,
                                f = '+')

  #Compute the average of M S_g M' for all g <= t_min
  avg_MSM_list <- purrr::map(.x = S_g_list[g_geq_gMin_index],
                             .f = function(S){ M %*% S %*% base::t(M) } )
  avg_MSM <- base::Reduce(x = avg_MSM_list,
                          f='+') / length(avg_MSM_list)

  #calculate the adjustment factor for the conservative variance
  N <- purrr::reduce(N_g_list,
                     sum)
  adjustmentFactor <- 1/N * base::t(betahat_g_sum) %*% avg_MSM %*% betahat_g_sum



  var_conservative <- seConservative^2
  if(var_conservative - adjustmentFactor < 0){
    warning("var_conservative is less than adjustmentFactor")
    se_adjusted <- 0
  }else{
    se_adjusted <- sqrt( var_conservative - adjustmentFactor )
  }

  if(!return_beta_sum){
    return(se_adjusted)
  }else{
    resultsList <- list(se = se_adjusted,
                        betahat_g_sum = betahat_g_sum,
                        avg_MSM = avg_MSM,
                        N= N)
    return(resultsList)
  }
}


#' @title create_A0_list
#' @description \code{create_A0_list} creates the list of A_0 matrices for Xhat corresponding with all possible
#' comparisons of cohorts before they are treated
#' @param g_list g_list
#' @param t_list t_list
#' @return A0_list list of A_0 matrices for Xhat corresponding with all possible comparisons of cohorts before they are treated

create_A0_list <- function(g_list,
                           t_list){

  createAtilde0_g <- function(g_index){
    g <- g_list[g_index]
    t_less_than_g_index <- which(t_list < g)

    #If no periods less than g, there are no restrictions, return null
    if(length(t_less_than_g_index) == 0){
      return(NULL)
    }

    #Otherwise, Atilde is a matrix with rows= # periods before g and T columns.
    # It has 1s on the diagonal and 0s elsewhere. I.e. Atilde selects rows of Y corresponding with t<g
    Atilde0_g <- matrix(0,
                        nrow = length(t_less_than_g_index),
                        ncol = length(t_list))
    diag(Atilde0_g) <- 1

    return(Atilde0_g)
  }

  Atilde0_list <- purrr::map(.x = 1:(length(g_list)-1),
                             .f = createAtilde0_g)

  A0_list <- purrr::map(.x = 1:length(Atilde0_list),
                        .f = function(i){
                          #A0_g is a block matrix that stacks the Atilde_g's (for g < gmax) and multiplies by 0 for all g' \neq g
                          Atilde_times_indicator_list <- purrr::map(.x = 1:length(Atilde0_list),
                                                                    .f = ~(i == .x) * Atilde0_list[[.x]]  )
                          A0 <- purrr::reduce(Atilde_times_indicator_list,
                                              rbind)
                          return(A0)
                        }  )
  if(length(A0_list) == 1){
    A0_gmax <- -A0_list[[1]]}
  else{
    A0_gmax <- -base::Reduce(x = A0_list,
                             f = '+')
  }

  A0_list[[length(A0_list) + 1]] <- A0_gmax
  return(A0_list)
}






sum_of_lists <- function(.l){
  #Given a list of lists all of the same length, this returns a list where the jth entry is the sum of the jth entry of all the lists
  if(length(.l) == 1){return(.l[[1]])}
  results <- purrr::reduce(.x = .l,
                           .f = function(l1,
                                         l2){
                             purrr::map2(l1,
                                         l2,
                                         .f = ~.x + .y)
                           }
  )
  return(results)
}

scalar_product_lists <- function(c,.l){
  #Takes a constant c and list .l, returns a list with c times each element of l
  results <- purrr::map(.x = .l,
                        .f = ~c * .x)
  return(results)
}

left_product_lists <- function(c,.l){
  #Takes a constant vector or matrix c and a list .l of conformable elements,
  # and returns a list with c %*% x for each element x of l
  results <- purrr::map(.x = .l,
                        .f = ~c %*% .x)
  return(results)
}


right_product_lists <- function(c,.l){
  #Takes a constant vector or matrix c and a list .l of conformable elements,
  # and returns a list with x %*% c for each element x of l
  results <- purrr::map(.x = .l,
                        .f = ~c %*% .x)
  return(results)
}


stack_rows_of_lists <- function(.l){
  #Given a list of lists all of the same length, where each inner list is a vector, this returns a list where the jth element is a matrix stacking all the jth elemtns
  if(length(.l) == 1){return(.l)}
  results <- purrr::reduce(.x = .l,
                           .f = function(l1,
                                         l2){
                             purrr::map2(l1,
                                         l2,
                                         .f = ~rbind(.x , .y))
                           }
  )
  return(results)
}

create_Atheta_list_for_ATE_tg <- function(t,
                                          g,
                                          g_list,
                                          t_list,
                                          N_g_list,
                                          use_last_treated_only = FALSE,
                                          showWarnings = TRUE){


  numPeriods <- length(t_list)
  #if(t < g & showWarnings){warning("t is less than g. ATE(t,g) is zero by assumption")}
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
    control_weights <- matrix(0, nrow =1, ncol = numPeriods)

    #If control_g is not a valid control, return 0
    if(! (control_g_index %in% control_cohort_indices[[treated_g_index]]) ){
      return(control_weights)
    }

    g_treated <- g_list[treated_g_index]
    N_g_control <- N_g_list[[control_g_index]]
    N_g_treated <- N_g_list[[treated_g_index]]

    t_control_index <- which(t_list == t)

    control_weights[t_control_index] <- - N_g_control / N_control_cohorts[treated_g_index]
    return(control_weights)
  }

  # createControlWeights <- function(control_g_index){
  #   #Create the control weights for control_g_index using all treated cohorts, then sum them
  #   controlWeights_g <- Reduce( x=  purrr::map(.x = eligible_cohort_index,
  #                                       ~createControlWeights_helper(control_g_index, .x)  ) ,
  #                               f= '+' )
  #   return(controlWeights_g)
  # }

  A_theta_treated_list <- purrr::map(.x = 1:length(g_list),
                                     A_theta_g_treated_fn)
  A_theta_control_list <- purrr::map(.x = 1:length(g_list),
                                     createControlWeights_helper )
  A_theta_list <- purrr::map2(.x = A_theta_treated_list,
                              .y = A_theta_control_list,
                              .f = ~ .x + .y)
  return(A_theta_list)
}



create_Atheta_list_for_event_study <- function(eventTime,
                                               g_list,
                                               t_list,
                                               N_g_list,
                                               use_last_treated_only = FALSE){

  #Create A_thetas for an ``event-study'' coefficient at lag eventTime
  # This is the average treatment effects for units eventTime periods after first being treated
  # The average is taken over all cohorts g such that there is some untreated cohort at g+eventTime
  # Cohorts are weighted by the cohort size (N_g)

  maxG <- max(g_list)
  eligible_cohort_index <- which( ((g_list + eventTime) < maxG ) & ((g_list + eventTime) <= max(t_list) ) )

  if(length(eligible_cohort_index) == 0){
    stop("There are no comparison cohorts for the given eventTime")
  }

  N_eligible <- base::Reduce(x = N_g_list[eligible_cohort_index],
                             sum)

  A_theta_lists <- purrr::map(.x = eligible_cohort_index,
                              .f = ~ scalar_product_lists(N_g_list[[.x]]/N_eligible ,
                                                          create_Atheta_list_for_ATE_tg(t =g_list[.x]+eventTime,
                                                                                        g = g_list[.x],
                                                                                        g_list = g_list,
                                                                                        t_list = t_list,
                                                                                        N_g_list = N_g_list,
                                                                                        use_last_treated_only = use_last_treated_only
                                                          )
                              )
  )

  if(length(eligible_cohort_index) == 1){
    A_theta_list <- A_theta_lists[[1]]
  }else{
    A_theta_list <- sum_of_lists(A_theta_lists)
  }
  return(A_theta_list)
}




create_Atheta_list_for_ATE_calendar_t <- function(t,
                                                  g_list,
                                                  t_list,
                                                  N_g_list,
                                                  use_last_treated_only = FALSE){

  treated_by_t_indices <- which(g_list <= t)
  N_total_treated <- sum( unlist(N_g_list[treated_by_t_indices]) )

  A_theta_lists <- purrr::map(.x = treated_by_t_indices,
                              .f = ~ scalar_product_lists(N_g_list[[.x]]/N_total_treated ,
                                                          create_Atheta_list_for_ATE_tg(t =t,
                                                                                        g = g_list[.x],
                                                                                        g_list = g_list ,
                                                                                        t_list = t_list,
                                                                                        N_g_list = N_g_list,
                                                                                        use_last_treated_only = use_last_treated_only
                                                          )
                              )
  )

  if(length(treated_by_t_indices) == 1){
    A_theta_list <- A_theta_lists[[1]]
  }else{
    A_theta_list <- sum_of_lists(A_theta_lists)
  }
  return(A_theta_list)
}




create_Atheta_list_for_ATE_cohort_g <- function(g,
                                                g_list,
                                                t_list,
                                                N_g_list,
                                                use_last_treated_only = FALSE){

  treated_period_indices <- which(t_list >= g & t_list < max(g_list))
  T_treated <- length( t_list[treated_period_indices] )

  A_theta_lists <- purrr::map(.x = treated_period_indices,
                              .f = ~ scalar_product_lists(1/T_treated ,
                                                          create_Atheta_list_for_ATE_tg(t =t_list[.x],
                                                                                        g = g,
                                                                                        g_list = g_list ,
                                                                                        t_list = t_list,
                                                                                        N_g_list = N_g_list,
                                                                                        use_last_treated_only = use_last_treated_only
                                                          )
                              )
  )

  if(T_treated == 1){
    A_theta_list <- A_theta_lists[[1]]
  }else{
    A_theta_list <- sum_of_lists(A_theta_lists)
  }

  return(A_theta_list)
}




create_Atheta_list_for_cohort_average_ATE <- function(g_list,
                                                      t_list,
                                                      N_g_list,
                                                      use_last_treated_only = FALSE){

  g_eligible_index <-  which((g_list < max(g_list)) & (g_list <= max(t_list)))

  N_total_eligible <- sum(unlist(N_g_list[g_eligible_index]))

  A_theta_lists <- purrr::map(.x = g_eligible_index,
                              .f = ~ scalar_product_lists(as.numeric(N_g_list[[.x]])/N_total_eligible,
                                                          create_Atheta_list_for_ATE_cohort_g(g = g_list[.x],
                                                                                              g_list = g_list ,
                                                                                              t_list = t_list,
                                                                                              N_g_list = N_g_list,
                                                                                              use_last_treated_only
                                                          )
                              )
  )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}




create_Atheta_list_for_calendar_average_ATE <- function(g_list,
                                                        t_list,
                                                        N_g_list,
                                                        use_last_treated_only = FALSE){

  t_eligible_index <-  which((t_list >= min(g_list)) & (t_list < max(g_list)))

  T_eligible <- length(t_eligible_index)

  A_theta_lists <- purrr::map(.x = t_eligible_index,
                              .f = ~ scalar_product_lists(1/T_eligible,
                                                          create_Atheta_list_for_ATE_calendar_t(t = t_list[.x],
                                                                                                g_list = g_list ,
                                                                                                t_list = t_list,
                                                                                                N_g_list = N_g_list,
                                                                                                use_last_treated_only = use_last_treated_only
                                                          )
                              )
  )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}


create_Atheta_list_for_simple_average_ATE <- function(g_list,
                                                      t_list,
                                                      N_g_list,
                                                      use_last_treated_only = FALSE){

  # Avoif Notes on CRAN
  g <- NULL

  #Create a df with all the (g,t) pairs for which ATE is identified
  gt_df <- purrr::cross_df( list(g = g_list, t = t_list) )
  gt_df <- gt_df %>% dplyr::filter(t >= g, t< max(g_list))

  #Join in N_g for each of these pairs
  gt_df <- dplyr::left_join( gt_df, data.frame(g= g_list, N_g = unlist(N_g_list)), by = "g")

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- sum(gt_df$N_g)

  A_theta_lists <- purrr::map(.x = 1:NROW(gt_df),
                              .f = ~ scalar_product_lists(gt_df$N_g[.x]/N_total,
                                                          create_Atheta_list_for_ATE_tg(t = gt_df$t[.x],
                                                                                        g= gt_df$g[.x],
                                                                                        g_list = g_list ,
                                                                                        t_list = t_list,
                                                                                        N_g_list = N_g_list,
                                                                                        use_last_treated_only = use_last_treated_only
                                                          )
                              )
  )


  A_theta_list <- sum_of_lists(A_theta_lists)

  return(A_theta_list)
}

#source(here("Code/refine-variance-estimates.R"))


calculate_full_vcv <- function(eventPlotResultsList, resultsDF){
  #Calculate A_theta_list - A_0_list %*% beta for each list
  #The reuslting list combined_A_list is a list of matrices of length |G| so that the vector of thetas is \sum combined_A_list[g] %*% Ybar
  combine_A_lists <- function(A_theta_list,
                              A_0_list,
                              beta){

    #Compute beta' %*% A_0_list
    A_0_beta_list <- left_product_lists(c = base::t(beta), .l = A_0_list)
    #Compute A_theta_list - A_0_list %*% beta
    combined_list <- sum_of_lists(list(A_theta_list, scalar_product_lists(-1,A_0_beta_list)) )
    return(combined_list)
  }
  combined_A_list <- purrr::map(.x = eventPlotResultsList,
                                .f = ~combine_A_lists(A_theta_list = .x$A_theta_list, A_0_list = .$A_0_list, beta = .$beta)  )

  combined_A_list <- stack_rows_of_lists(combined_A_list)

  #The newman vcv is \sum 1/N_g * combined_A * S * combined_A'
  vcv_neyman_terms_list <-
    purrr::pmap(.l = list(S_g = eventPlotResultsList[[1]]$S_g_list, A = combined_A_list, N_g = eventPlotResultsList[[1]]$N_g_list),
                .f = function(S_g, A, N_g){ return((1/N_g)* A %*% S_g %*% base::t(A) ) } )
  vcv_neyman <- base::Reduce(f = '+', x = vcv_neyman_terms_list)


  stacked_betahat_g_sum <- base::Reduce(x = purrr::map(.x = eventPlotResultsList, .f = ~base::t(.x$betahat_g_sum)),
                                        f = rbind)

  vcv_adjustment <- 1/eventPlotResultsList[[1]]$N * stacked_betahat_g_sum %*% eventPlotResultsList[[1]]$avg_MSM %*% base::t(stacked_betahat_g_sum)

  vcv <- vcv_neyman - vcv_adjustment

  return(list(vcv = vcv, vcv_neyman = vcv_neyman))
}



processDF <- function(df, i, g, t, y){

  #This function processes the df inputted to staggered (or staggered_cs/sa)
  # It checks that the columns in the user-inputted values of  i,g,t,y are actually in the data
  # It also renames these columns to "i", "g", "t", "y"

  # Let's make sure we have columns with name i, t, y and g
  colnames_df <- colnames(df)
  if(!i %in% colnames_df){
    stop(paste0("There is no column ", i, " in the data. Thus, we are not able to find the unit identifier variable."))
  }
  if(!t %in% colnames_df){
    stop(paste0("There is no column ", t, " in the data. Thus, we are not able to find the time identifier variable."))
  }
  if(!g %in% colnames_df){
    stop(paste0("There is no column ", g, " in the data. Thus, we are not able to find the group identifier variable."))
  }
  if(!y %in% colnames_df){
    stop(paste0("There is no column ", y, " in the data. Thus, we are not able to find the outcome variable."))
  }

  # Sanity checks
  if(i %in% c("g", "t", "y" )){
    stop(paste0("Unit identifier cannot be labeled g, t, or y"))
  }

  if(t %in% c("i","y", "g" )){
    stop(paste0("Time identifier cannot be labeled i, g, or y"))
  }

  if(g %in% c("i", "t" ,"y" )){
    stop(paste0("Group identifier cannot be labeled i, t, or y"))
  }


  # Re-label i, t, g, y
  if(i != "i"){
    df[,"i"] <- df[,i]
  }

  if(t != "t"){
    df[, "t"] <- df[,t]
  }

  if(g != "g"){
    df[, "g"] <-  df[,g]
  }

  if(y != "y"){
    df[, "y"] <- df[,y]
  }
  return(df)
}


#' @useDynLib staggered
#' @importFrom magrittr "%>%"
#' @import Rcpp
#' @title Calculate the efficient adjusted estimator in staggered rollout designs
#' @description This functions calculates the efficient estimator for staggered rollout designs proposed by Roth and Sant'Anna.
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param i The name of column containing the individual (cross-sectional unit) identifier. Default is "i".
#' @param t The name of the column containing the time periods. Default is "t".
#' @param g The name of the column containing the first period when a particular observation is treated, with Inf denoting never treated. Default is "g".
#' @param y The name of the column containing the outcome variable. Default is "y".
#' @param estimand The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. "eventstudy" returns the average effect at the ''event-time'' given in the parameter EventTime.  The parameter can be left blank if a custom parameter is provided in A_theta_list. The argument is not case-sensitive.
#' @param A_theta_list This parameter allows for specifying a custom estimand, and should be left as NULL if estimand is specified. It is a list of matrices A_theta_g so that the parameter of interest is sum_g A_theta_g Ybar_g, where Ybar_g = 1/N sum_i Y_i(g)
#' @param A_0_list This parameter allow for specifying the matrices used to construct the Xhat vector of pre-treatment differences. If left NULL, the default is to use the scalar set of controls used in Callaway and Sant'Anna. If use_DiD_A0 = FALSE, then it uses the full vector possible comparisons of (g,g') in periods t<g,g'.
#' @param eventTime If using estimand = "eventstudy", specify what eventTime you want the event-study parameter for. The default is 0, the period in which treatment occurs. If a vector is provided, estimates are returned for all the event-times in the vector.
#' @param beta A coefficient to use for covariate adjustment. If not specified, the plug-in optimal coefficient is used. beta =0 corresponds with the simple difference-in-means. beta = 1 corresponds with the Callaway and Sant'Anna estimator when using the default value of use_DiD_A0 = TRUE.
#' @param use_DiD_A0 If this parameter is true, then Xhat corresponds with the scalar used by Callaway and Sant'Anna, so the Callaway and Sant'Anna estimator corresponds with beta=1. If it is false, the Xhat is a vector with all possible comparisons of pairs of cohorts before either is treated. The latter option should only be used when the number of possible comparisons is small relative to sample size.
#' @param return_full_vcv If this is true and estimand = "eventstudy", then the function returns a list containing the full variance-covariance matrix for the event-plot estimates in addition to the usual dataframe with the estimates
#' @param return_matrix_list If true, the function returns a list of the A_0_list and A_theta_list matrices along with betastar. This is used for internal recursive calls to calculate the variance-covariance matrix, and will generally not be needed by the end-user. Default is False.
#' @param use_last_treated_only If true, then A_0_list and A_theta_list are created to only make comparisons with the last treated cohorts (as suggested by Sun and Abraham), rather than using not-yet-treated units as comparisons. If set to TRUE (and use_DiD_A0 = TRUE), then beta=1 corresponds with the Sun and Abraham estimator.
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.
#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman). (If return_matrix_list = TRUE, it likewise returns a list containing lists of matrices used in the vcv calculation.)
#' @references
#' \cite{Roth, Jonatahan, and Sant'Anna, Pedro H. C. (2021),
#'   'Efficient Estimation for Staggered Rollout Designs', arXiv: 2102.01291, \url{https://arxiv.org/abs/2102.01291}.}
#' @examples
#'
#' # Load some libraries
#' library(dplyr)
#' library(purrr)
#' library(MASS)
#' set.seed(1234)
#' # load the officer data and subset it
#' df <- pj_officer_level_balanced
#' group_random <- sample(unique(df$assigned), 3)
#' df <- df[df$assigned %in% group_random,]
#' # Calculate efficient estimator for the simple weighted average
#' staggered(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "simple")
#' # Calculate efficient estimator for the cohort weighted average
#' staggered(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "cohort")
#' # Calculate efficient estimator for the calendar weighted average
#' staggered(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "calendar")
#' # Calculate event-study coefficients for the first 24 months
#' # (month 0 is instantaneous effect)
#' eventPlotResults <- staggered(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "eventstudy",
#'   eventTime = 0:23)
#' eventPlotResults %>% head()
#'
#' @export
staggered <- function(df,
                      i = "i",
                      t = "t",
                      g = "g",
                      y = "y",
                      estimand = NULL,
                      A_theta_list = NULL,
                      A_0_list = NULL,
                      eventTime = 0,
                      beta = NULL,
                      use_DiD_A0 = ifelse(is.null(A_0_list),
                                          TRUE,
                                          FALSE),
                      return_full_vcv = FALSE,
                      return_matrix_list = FALSE,
                      use_last_treated_only = FALSE,
					            compute_fisher = FALSE,
					            num_fisher_permutations = 500,
					            skip_data_check = FALSE){

  #Process the inputted df by checking the inputted columns and renaming to i,t,g,y, and balancing on (i,t)
  #We skip this if skip_data_check = TRUE
  if(!skip_data_check){
    df <- processDF(df,
                    i=i,
                    g=g,
                    t=t,
                    y=y)
    #Balance the panel (and throw a warning if original panel is unbalanced)
    df <- balance_df(df = df)
  }


  #  Compute number of units per cohort
  cohort_size <- base::table(df$g)/base::length(base::table(df$t))
  # Flag for singleton cohorts
  flag_singleton <- as.numeric(names(cohort_size[(cohort_size==1)]))
  # Drop cohorts which are singleton
  l_flag1 <- base::length(flag_singleton)
  if(l_flag1 > 0){
    gpaste <-  paste(flag_singleton, collapse=", ")
    if(l_flag1==1){
      base::warning(paste0("The treatment cohort g = ", gpaste, " has a single cross-sectional unit. We drop this cohort."))
    } else {
      base::warning(paste0("The treatment cohorts g = ", gpaste, " have a single cross-sectional unit only. We drop these cohorts."))
    }

   df <- df[(df$g %in% flag_singleton) == FALSE,]
  }


  #  If estimand is provided, force to be lower-case (allowing for non-case sensitive inputs)
  if(!is.null(estimand)){
    estimand <- tolower(estimand)
  }

  #If eventTime is a vector, call staggered for each event-time and combine the results
  #Add the variable eventTime to the data frame
  if(length(eventTime) > 1){

    if(estimand != "eventstudy"){
      stop("You provided a vector fpr eventTime but estimand is not set to 'eventstudy'. Did you mean to set estimand = 'eventstudy'?")
    }

    eventPlotResultsList <-
      purrr::map(.x = eventTime,
                 .f = ~staggered(df = df,
                                 estimand = estimand,
                                 A_theta_list = A_theta_list,
                                 A_0_list = A_0_list,
                                 eventTime = .x,
                                 beta = beta,
                                 use_DiD_A0 = use_DiD_A0,
                                 return_matrix_list = TRUE,
                                 use_last_treated_only = use_last_treated_only,
                                 compute_fisher = compute_fisher,
                                 skip_data_check = T))

    resultsDF <- purrr::reduce(.x = purrr::map(.x = eventPlotResultsList, .f = ~ .x$resultsDF),
                               .f = dplyr::bind_rows)

    #Add in eventTimes
    resultsDF$eventTime <- eventTime

    if(return_full_vcv){

      vcvs <- calculate_full_vcv(eventPlotResultsList = eventPlotResultsList,
                                 resultsDF = resultsDF)



      resultsList <- list(resultsDF = resultsDF, vcv = vcvs$vcv, vcv_neyman = vcvs$vcv_neyman)

      #Create stacked beta for the
      return(resultsList)
    }else{
      return(resultsDF)
    }
  }
  g_level_summaries <- compute_g_level_summaries(df, is_balanced = TRUE)
  Ybar_g_list <- g_level_summaries$Ybar_g_List
  S_g_list <- g_level_summaries$S_g_List
  N_g_list <- g_level_summaries$N_g_List
  g_list <- g_level_summaries$g_list
  t_list <- g_level_summaries$t_list


  #If estimand is provided, calculate the appropriate A_theta_list
  if(!is.null(estimand)){
    if(estimand == "simple"){
      A_theta_list <- create_Atheta_list_for_simple_average_ATE(g_list = g_list,
                                                                t_list = t_list,
                                                                N_g_list = N_g_list,
                                                                use_last_treated_only = use_last_treated_only)
    }else if(estimand == "cohort"){
      A_theta_list <- create_Atheta_list_for_cohort_average_ATE(g_list = g_list,
                                                                t_list = t_list,
                                                                N_g_list = N_g_list,
                                                                use_last_treated_only = use_last_treated_only)
    }else if(estimand == "calendar"){
      A_theta_list <- create_Atheta_list_for_calendar_average_ATE(g_list = g_list,
                                                                  t_list = t_list,
                                                                  N_g_list = N_g_list,
                                                                  use_last_treated_only = use_last_treated_only)
    }else if(estimand == "eventstudy"){
      A_theta_list <- create_Atheta_list_for_event_study(eventTime = eventTime,
                                                         g_list = g_list,
                                                         t_list = t_list,
                                                         N_g_list = N_g_list,
                                                         use_last_treated_only = use_last_treated_only)
    }
  }
  #If no valid estimand is provided and no A_theta_list, throw and error
  if(is.null(A_theta_list)){
    stop("Estimand must be one of simple, cohort, calendar, or eventstudy; or custom A_theta_list must be provided")
  }

  #Create A_0_list if a custom A_0_list is not provided
  if(is.null(A_0_list) & (use_DiD_A0==FALSE)){
    A_0_list <- create_A0_list(g_list = g_list,
                               t_list = t_list)
  }

  #If use_DiD_A0, use only the A0's associated with the DiD estimand
  if(use_DiD_A0){

    if(is.null(estimand)){
      stop("If use_DiD_A0 = TRUE, you must provide an estimand.")
    }

    if(estimand == "simple"){
      A_0_list <- create_A0_list_for_simple_average_ATE(g_list = g_list,
                                                        t_list = t_list,
                                                        N_g_list = N_g_list,
                                                        use_last_treated_only = use_last_treated_only)
    }else if(estimand == "cohort"){
      A_0_list <- create_A0_list_for_cohort_average_ATE(g_list = g_list,
                                                        t_list = t_list,
                                                        N_g_list = N_g_list,
                                                        use_last_treated_only = use_last_treated_only)
    }else if(estimand == "calendar"){
      A_0_list <- create_A0_list_for_calendar_average_ATE(g_list = g_list,
                                                          t_list = t_list,
                                                          N_g_list = N_g_list,
                                                          use_last_treated_only = use_last_treated_only)
    }else if(estimand == "eventstudy"){
      A_0_list <- create_A0_list_for_event_study(eventTime = eventTime,
                                                 g_list = g_list,
                                                 t_list = t_list,
                                                 N_g_list = N_g_list,
                                                 use_last_treated_only = use_last_treated_only)
    }
  }


  Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) ,
                           .f = function(A0,S,N){
                             return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , base::t(A0) ) )
                           }
  )
  #  Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list) , .f = function(A0,S,N){ return(1/N * A0 %*%S %*% base::t(A0) )  } )

  #Save the user-inputted beta (used in FRT call)
  user_input_beta <- beta

  if(is.null(beta)){
    beta <- compute_Betastar(Ybar_g_list,
                             A_theta_list,
                             A_0_list,
                             S_g_list,
                             N_g_list,
                             Xvar_list = Xvar_list)
  }

  #If beta =0, convert beta to the appropriate length
  if(length(beta) == 1 && beta == 0){
    beta = matrix(0, dim(A_0_list[[1]])[1])
  }

  thetahat <- compute_Thetahat_beta(beta = beta,
                                    Ybar_g_list,
                                    A_theta_list,
                                    A_0_list,
                                    S_g_list,
                                    N_g_list,
                                    Xvar_list = Xvar_list)
  se_neyman <- compute_se_Thetahat_beta_conservative(beta = beta,
                                                     Ybar_g_list,
                                                     A_theta_list,
                                                     A_0_list,
                                                     S_g_list,
                                                     N_g_list,
                                                     Xvar_list = Xvar_list)
  seResults <- compute_se_Thetahat_beta(beta = beta,
                                        Ybar_g_list,
                                        A_theta_list,
                                        A_0_list,
                                        S_g_list,
                                        N_g_list,
                                        g_list,
                                        t_list,
                                        Xvar_list = Xvar_list,
                                        return_beta_sum = TRUE
  )

  se <- seResults$se

  resultsDF <- data.frame(estimate = thetahat,
                          se = se,
                          se_neyman = se_neyman)

  ## Do FRT, if specified
  permuteTreatment <- function(df,i_g_table, seed){
    #This function takes a data.frame with columns i and g, and permutes the values of g assigned to i
    # The input i_g_table has the unique combinations of (i,g) in df, and is calculated outside for speed improvements

    #Draw a random permutation of the elements of first_period_df
    set.seed(seed)
    n = NROW(i_g_table)
    randIndex <-
      sample.int(n = n,
                 size = n,
                 replace = F)

    #Replace first_period_df$g with a permuted version based on randIndex
    i_g_table$g <- i_g_table$g[randIndex]

    #Merge the new treatment assignments back with the original
    df$g <- NULL
    df <- dplyr::left_join(df,
                           i_g_table,
                           by = c("i"))

    return(df)
  }

  if(compute_fisher){

    #Find unique pairs of (i,g). This will be used for computing the permutations
    # i_g_table <- df %>%
    #              dplyr::filter(t == min(t)) %>%
    #              dplyr::select(i,g)

    i_g_table <- df %>%
      dplyr::filter(t == min(t))
    i_g_table <- i_g_table[,c("i","g")]

    #Now, we compute the FRT for each seed, permuting treatment for each one
      #We catch any errors in the FRT simulations, and throw a warning if at least one has an error (using the remaining draws to calculate frt)
    FRTResults <-
      purrr::map(.x = 1:num_fisher_permutations,
                 .f = purrr::possibly(
                   .f =~ staggered::staggered(df = permuteTreatment(df, i_g_table, seed = .x),
                                              estimand = NULL,
                                              beta = user_input_beta,
                                              A_theta_list = A_theta_list,
                                              A_0_list = A_0_list,
                                              eventTime = eventTime,
                                              return_full_vcv = F,
                                              return_matrix_list = F,
                                              compute_fisher = F,
                                              skip_data_check = T) %>% mutate(seed = .x),
                   otherwise = NULL)
      ) %>%
      purrr::discard(base::is.null) %>%
      purrr::reduce(.f = dplyr::bind_rows)

    successful_frt_draws <- NROW(FRTResults)
    if(successful_frt_draws < num_fisher_permutations){
      warning("There was an error in at least one of the FRT simulations. Removing the problematic draws.")
    }

    resultsDF$fisher_pval <- mean( abs(resultsDF$estimate/resultsDF$se) < abs(FRTResults$estimate/FRTResults$se) )
    resultsDF$fisher_pval_se_neyman <- mean( abs(resultsDF$estimate/resultsDF$se_neyman) < abs(FRTResults$estimate/FRTResults$se_neyman) )
    resultsDF$num_fisher_permutations <- successful_frt_draws

  }
  if(!return_matrix_list){
    #If return_matrix_list is not specified, then we return results DF unless return_full_vcv = TRUE
    if(!return_full_vcv){
      return(resultsDF)
    }else{
      resultsList <- list(resultsDF = resultsDF,
                          vcv = as.matrix(se^2),
                          vcv_neyman = as.matrix(se_neyman^2))

      return(resultsList)
    }


  }else{
    resultsList <- list(resultsDF = resultsDF,
                        A_theta_list = A_theta_list,
                        A_0_list = A_0_list,
                        beta = beta,
                        betahat_g_sum = seResults$betahat_g_sum,
                        avg_MSM = seResults$avg_MSM,
                        S_g_list = S_g_list,
                        N_g_list = N_g_list,
                        N = seResults$N)

    return(resultsList)
  }


}




#' @title Calculate the Callaway & Sant'Anna (2020) estimator for staggered rollouts
#' @description This functions calculates the Callaway & Sant'Anna (2020) estimator for staggered rollout designs using not-yet-treated units (including never-treated, if available) as controls.
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param i The name of column containing the individual (cross-sectional unit) identifier. Default is "i".
#' @param t The name of the column containing the time periods. Default is "t".
#' @param g The name of the column containing the first period when a particular observation is treated, with Inf denoting never treated. Default is "g".
#' @param y The name of the column containing the outcome variable. Default is "y".
#' @param estimand The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. "eventstudy" returns the average effect at the ''event-time'' given in the parameter EventTime.  The parameter can be left blank if a custom parameter is provided in A_theta_list. The argument is not case-sensitive.
#' @param A_theta_list This parameter allows for specifying a custom estimand, and should be left as NULL if estimand is specified. It is a list of matrices A_theta_g so that the parameter of interest is sum_g A_theta_g Ybar_g, where Ybar_g = 1/N sum_i Y_i(g)
#' @param A_0_list This parameter allow for specifying the matrices used to construct the Xhat vector of pre-treatment differences. If left NULL, the default is to use the scalar set of controls used in Callaway and Sant'Anna. If use_DiD_A0 = FALSE, then it uses the full vector possible comparisons of (g,g') in periods t<g,g'.
#' @param eventTime If using estimand = "eventstudy", specify what eventTime you want the event-study parameter for. The default is 0, the period in which treatment occurs. If a vector is provided, estimates are returned for all the event-times in the vector.
#' @param return_full_vcv If this is true and estimand = "eventstudy", then the function returns a list containing the full variance-covariance matrix for the event-plot estimates in addition to the usual dataframe with the estimates
#' @param return_matrix_list If true, the function returns a list of the A_0_list and A_theta_list matrices along with betastar. This is used for internal recursive calls to calculate the variance-covariance matrix, and will generally not be needed by the end-user. Default is False.
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.

#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman). (If return_matrix_list = TRUE, it likewise returns a list containing lists of matrices used in the vcv calculation.)
#' @references
#'   \cite{Callaway, Brantly, and Sant'Anna, Pedro H. C. (2020),
#'   'Difference-in-Differences with Multiple Time Periods', Forthcoming at the Journal of Econometrics,
#'   \doi{10.1016/j.jeconom.2020.12.001}.}
#' @examples
#' # Load some libraries
#' library(dplyr)
#' library(purrr)
#' library(MASS)
#' set.seed(1234)
#' # load the officer data and subset it
#' df <- pj_officer_level_balanced
#' group_random <- sample(unique(df$assigned), 3)
#' df <- df[df$assigned %in% group_random,]
#' # We modify the data so that the time dimension is named t,
#' # the period of treatment is named g,
#' # the outcome is named y,
#' # and the individual identifiers are named i
#'# (this allow us to use default arguments on \code{staggered_cs}).
#' df <- df %>% rename(t = period, y = complaints, g = first_trained, i = uid)
#' # Calculate Callaway and Sant'Anna estimator for the simple weighted average
#' staggered_cs(df = df, estimand = "simple")
#' # Calculate Callaway and Sant'Anna estimator for the cohort weighted average
#' staggered_cs(df = df, estimand = "cohort")
#' # Calculate Callaway and Sant'Anna estimator for the calendar weighted average
#' staggered_cs(df = df, estimand = "calendar")
#' # Calculate Callaway and Sant'Anna event-study coefficients for the first 24 months
#' # (month 0 is instantaneous effect)
#' eventPlotResults <- staggered_cs(df = df, estimand = "eventstudy", eventTime = 0:23)
#' eventPlotResults %>% head()
#'
#' @export
staggered_cs <- function(df,
                         i = "i",
                         t = "t",
                         g = "g",
                         y = "y",
                         estimand = NULL,
                         A_theta_list = NULL,
                         A_0_list = NULL,
                         eventTime = 0,
                         return_full_vcv = FALSE,
                         return_matrix_list = FALSE,
                         compute_fisher = FALSE,
                         num_fisher_permutations = 500,
                         skip_data_check = FALSE){

  if(!skip_data_check){
    df <- processDF(df,
                    i=i,
                    g=g,
                    t=t,
                    y=y)
    #Balance the panel (and throw a warning if original panel is unbalanced)
    df <- balance_df(df = df)
  }


  #Drop units who has g= < min(t), since ATT(t,g) is not identified for these units
  TreatedBeforeMinT <- df$g <= min(df$t)
  if(sum(TreatedBeforeMinT) > 0 ){
    df <- df[!TreatedBeforeMinT, ]
    warning("Dropping units who were treated in the first period or earlier, since CS estimator is not defined (and ATT(t,g) not identified under parallel trends).")
  }

  results <- staggered(df = df,
                       estimand = estimand,
                       A_theta_list = A_theta_list,
                       eventTime = eventTime,
                       beta = 1,
                       use_DiD_A0 = TRUE,
                       use_last_treated_only = FALSE,
                       return_full_vcv = return_full_vcv,
                       return_matrix_list = return_matrix_list,
                       compute_fisher = compute_fisher,
                       num_fisher_permutations = num_fisher_permutations,
                       skip_data_check = skip_data_check)

  return(results)

}




#' @title Calculate the Sun & Abraham (2020) estimator for staggered rollouts
#' @description This functions calculates the Sun & Abraham (2020) estimator for staggered rollout designs using last-treated-treated units (never-treated, if availabe) as controls.
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param i The name of column containing the individual (cross-sectional unit) identifier. Default is "i".
#' @param t The name of the column containing the time periods. Default is "t".
#' @param g The name of the column containing the first period when a particular observation is treated, with Inf denoting never treated. Default is "g".
#' @param y The name of the column containing the outcome variable. Default is "y".
#' @param estimand The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. "eventstudy" returns the average effect at the ''event-time'' given in the parameter EventTime.  The parameter can be left blank if a custom parameter is provided in A_theta_list. The argument is not case-sensitive.
#' @param A_theta_list This parameter allows for specifying a custom estimand, and should be left as NULL if estimand is specified. It is a list of matrices A_theta_g so that the parameter of interest is sum_g A_theta_g Ybar_g, where Ybar_g = 1/N sum_i Y_i(g)
#' @param A_0_list This parameter allow for specifying the matrices used to construct the Xhat vector of pre-treatment differences. If left NULL, the default is to use the scalar set of controls used in Callaway and Sant'Anna. If use_DiD_A0 = FALSE, then it uses the full vector possible comparisons of (g,g') in periods t<g,g'.
#' @param eventTime If using estimand = "eventstudy", specify what eventTime you want the event-study parameter for. The default is 0, the period in which treatment occurs. If a vector is provided, estimates are returned for all the event-times in the vector.
#' @param return_full_vcv If this is true and estimand = "eventstudy", then the function returns a list containing the full variance-covariance matrix for the event-plot estimates in addition to the usual dataframe with the estimates
#' @param return_matrix_list If true, the function returns a list of the A_0_list and A_theta_list matrices along with betastar. This is used for internal recursive calls to calculate the variance-covariance matrix, and will generally not be needed by the end-user. Default is False.
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.

#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman). (If return_matrix_list = TRUE, it likewise returns a list containing lists of matrices used in the vcv calculation.)
#' @references
#'   \cite{Sun, Liyang, and Abraham, Sarah (2020),
#'   'Estimating dynamic treatment effects in event studies with heterogeneous treatment effects', Forthcoming at the Journal of Econometrics,
#'   \doi{10.1016/j.jeconom.2020.09.006}.}
#' @examples

#' # Load some libraries
#' library(dplyr)
#' library(purrr)
#' library(MASS)
#' set.seed(1234)
#' # load the officer data and subset it
#' df <- pj_officer_level_balanced
#' group_random <- sample(unique(df$assigned), 3)
#' df <- df[df$assigned %in% group_random,]
#' # We modify the data so that the time dimension is named t,
#' # the period of treatment is named g,
#' # the outcome is named y,
#' # and the individual identifiers are named i
#' #  (this allow us to use default arguments on \code{staggered_cs}).
#' df <- df %>% rename(t = period, y = complaints, g = first_trained, i = uid)
#' # Calculate Sun and Abraham estimator for the simple weighted average
#' staggered_sa(df = df, estimand = "simple")
#' # Calculate Sun and Abraham estimator for the cohort weighted average
#' staggered_sa(df = df, estimand = "cohort")
#' # Calculate Sun and Abraham estimator for the calendar weighted average
#' staggered_sa(df = df, estimand = "calendar")
#' # Calculate Sun and Abraham event-study coefficients for the first 24 months
#' # (month 0 is instantaneous effect)
#' # eventPlotResults <- staggered_sa(df = df, estimand = "eventstudy", eventTime = 0:23)
#' # eventPlotResults %>% head()
#'
#' @export
staggered_sa <- function(df,
                         i = "i",
                         t = "t",
                         g = "g",
                         y = "y",
                         estimand = NULL,
                         A_theta_list = NULL,
                         A_0_list = NULL,
                         eventTime = 0,
                         return_full_vcv = FALSE,
                         return_matrix_list = FALSE,
                         compute_fisher = FALSE,
                         num_fisher_permutations = 500,
                         skip_data_check = FALSE){


  if(!skip_data_check){
    df <- processDF(df,
                    i=i,
                    g=g,
                    t=t,
                    y=y)
    #Balance the panel (and throw a warning if original panel is unbalanced)
    df <- balance_df(df = df)
  }

  #Drop units who has g= < min(t), since ATT(t,g) is not identified for these units
  TreatedBeforeMinT <- df$g <= min(df$t)
  if(sum(TreatedBeforeMinT) > 0 ){
    df <- df[!TreatedBeforeMinT, ]
    warning("Dropping units who were treated in the first period or earlier, since SA estimator is not defined (and ATT(t,g) not identified under parallel trends).")
  }

  results <- staggered(df = df,
                       estimand = estimand,
                       A_theta_list = A_theta_list,
                       eventTime = eventTime,
                       beta = 1,
                       use_DiD_A0 = TRUE,
                       use_last_treated_only = TRUE,
                       return_full_vcv = return_full_vcv,
                       return_matrix_list = return_matrix_list,
                       compute_fisher = compute_fisher,
                       num_fisher_permutations = num_fisher_permutations,
                       skip_data_check = skip_data_check)

  return(results)

}
