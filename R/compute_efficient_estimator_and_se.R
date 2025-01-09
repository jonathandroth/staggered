# ----------------------------------------------------------------------
# Misc helpers
# ----------------------------------------------------------------------

#' @title Calculate group level summary statistics
#' @description This function computes the mean-vector and covariance matrix of the outcomes for each cohort, where a cohort g is a group of units first treated in period g
#' @param df A data table containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param is_balanced If true, the df has previously been balanced so this does not need to be done internally.
#' @return Y_bar_list A list of the means of the outcomes for each cohort g
#' @return S_g_list A list of covariance matrices for the outcomes for each cohort g
#' @return N_g_DT A data table of the number of observations for each cohort g
#' @return g_list A list of when the cohorts were first treated
#' @return t_list A list of the the time periods for the outcome. The vector of outcomes corresponds with this order.
compute_g_level_summaries <- function(df, is_balanced = TRUE){
  N_g <- NULL

  #Balance the panel (and throw a warning if original panel is unbalanced)
  if(!is_balanced){
    df <- balance_df(df)
  }
  g_list <- base::sort(df[,base::unique(g)])
  t_list <- base::sort(df[,base::unique(t)])

  #Reshape so that Y_{it} is a column for all t
  #NB: Columns with data.table::dcast will be ordered
  #by unique values of t, per the docs.
  df <- data.table::dcast(df, i + g ~ t, value.var="y")[,!"i"]

  #Sort by g
  data.table::setkey(df, g)
  data.table::setindex(df, g)

  # Helper fun for covariance matrix
  compute_vcov_g <- function(group) {
      stats::var(df[.(group),!"g"])
  }

  Ybar_g_DT <- df[,base::lapply(.SD, base::mean), by="g"]
  N_g_DT    <- df[,.(N_g=.N),by="g"]
  S_g_List  <- base::lapply(g_list, compute_vcov_g)


  return(list(Ybar_g_List = base::asplit(Ybar_g_DT[,!"g"], 1),
              S_g_List    = S_g_List,
              N_g_DT      = N_g_DT,
              g_list      = g_list,
              t_list      = t_list))
}


balance_df <- function(df){
  numPeriods_i <- NULL

  # This function creates a balanced panel as needed for our analysis
  # -----------------------------------------------------------------
  #
  # - It first checks if rows of the data are uniquely characterized by (i,t);
  #   if there are multiple observations per (i,t), it throws an error
  # - It also removes observations with missing y
  # - Last, it removes observations i for which data is not available for all t

  # Sort by (i, t)
  data.table::setkey(df, i, t)

  # Check (i,t) is a unique identifier
  if(base::any(base::duplicated(df[,c("i","t")]))){
    base::stop("There are multiple observations with the same (i,t) values. The panel should have a unique outcome for each (i,t) value.")
  }

  # Check if there are missign values for y, and remove them if so
  if ( base::any(base::is.na(df$y)) ){
    df <- df[!base::is.na(df$y)]
  }

  # Check if panel is balanced. If not, drop the unbalanced observations and throw a warning
  df[,numPeriods_i := base::length(t), by=.(i)]
  numPeriods <- base::length(df[,base::unique(t)])
  if(df[,base::min(numPeriods_i)] == numPeriods){
    return(df[,!c('numPeriods_i')])
  }else{
    base::warning("Panel is unbalanced (or has missing values for some observations). Dropping observations with missing values of Y_{it} for some time periods. If you wish to include these observations, provide staggered with a df with imputed outcomes.")
    df <- df[numPeriods_i==numPeriods,!c('numPeriods_i')]
    return(df)
  }
}


compute_Thetahat0 <- function(Ybar_g_list,
                              A_theta_list){

  A_theta_Ybar_list <- purrr::map2(.x = Ybar_g_list,
                                   .y = A_theta_list,
                                   .f = ~.y %*% .x)
  Thetahat0 <- purrr::reduce(.x = A_theta_Ybar_list,
                             .f = base::sum)
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

  if(base::is.null(Xvar_list)){
    # Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, N_g_list),
    #                          .f = function(A0, S, N){ return(1/N * A0 %*% S %*% base::t(A0) ) } )
    Xvar_list <- purrr::pmap(.l = list(A_0_list,
                                       S_g_list, N_g_list),
                             .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , base::t(A0) ) ) }
    )
  }

  # Xvar <- purrr::reduce(.x = Xvar_list, .f = sum)
  Xvar <- base::Reduce(f = '+',
                       x= Xvar_list)

  X_theta_cov_list <- purrr::pmap(.l = list(A_0_list,
                                            A_theta_list,
                                            S_g_list,
                                            N_g_list),
                                  .f = function(A0,A_theta,S,N){ return(1/N * A0 %*% S %*% base::t(A_theta) ) }
  )

  # X_theta_cov <- purrr::reduce(.x = X_theta_cov_list, .f = sum)
  X_theta_cov <- base::Reduce(x = X_theta_cov_list, f = '+')

  # betastar <- solve(Xvar) %*% X_theta_cov
  # betastar <- MASS::ginv(Xvar) %*% X_theta_cov
  # betastar <- solve_least_squares_svd(Xvar,X_theta_cov)
  betastar <- solve_least_squares_normal(Xvar,X_theta_cov)  # fast method of solving (Xvar)^-1 X_theta_cov
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
    # Xvar_list <- pmap(.l = list(A_0_list, S_g_list, N_g_list) ,
    #                   .f = function(A0,S,N){ return(1/N * A0 %*% S %*% base::t(A0) ) } )
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
    base::warning("Calculated variance is less than 0. Setting to 0.")
    se_neyman <- 0
  }else{
    se_neyman <- base::sqrt(varhat_conservative)
  }
  return(se_neyman)
}


computeGMin <- function(A_theta_list,
                        g_list){

  A_theta_is_nonzero <- purrr::map_lgl(.x = A_theta_list,
                                       .f = ~base::max(base::abs(.x)) != 0 )

  min_nonzero_index <- base::min( base::which(A_theta_is_nonzero) )
  g_min <- g_list[min_nonzero_index]

  if(base::length(g_min) == 0){
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
                                     return_beta_sum = FALSE,
                                     gMin = NULL){

  # This function computes the standard error, using the version sigma_**
  # that adjust for pre-treatment covariates
  #
  # If return_beta_sum = TRUE, it returns a list wit
  seConservative <- compute_se_Thetahat_beta_conservative(beta,
                                                          Ybar_g_list,
                                                          A_theta_list,
                                                          A_0_list,
                                                          S_g_list,
                                                          N_g_list,
                                                          Xvar_list =  Xvar_list)
  if(base::is.null(gMin)){
    gMin <- computeGMin(A_theta_list = A_theta_list,
                        g_list = g_list)

  }

  g_geq_gMin_index <- base::which(g_list >= gMin) #which indices of g are >= gMin

  tMin <- base::min(t_list)
  tMax <- base::max(t_list)

  #Exit if first cohort is always treated
  if(gMin <= tMin){
    if(!return_beta_sum){
      return(seConservative)
    }else{
      resultsList <- list(se = seConservative,
                          betahat_g_sum = 0,
                          avg_MSM = 0,
                          N= purrr::reduce(N_g_list,
                                           base::sum))
      return(resultsList)
    }
  }

  #Create matrix M that selects the rows of S_g correspondign with t< g_min
  M <- matrix(0, nrow = gMin -tMin,
              ncol = base::NCOL(S_g_list[[1]]) )
  base::diag(M) <- 1

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
                             .f = ~(M %*% .x %*% base::t(M)) )
  avg_MSM <- base::Reduce(x = avg_MSM_list,
                          f='+') / base::length(avg_MSM_list)

  #calculate the adjustment factor for the conservative variance
  N <- purrr::reduce(N_g_list,
                     base::sum)
  adjustmentFactor <- 1/N * base::t(betahat_g_sum) %*% avg_MSM %*% betahat_g_sum

  var_conservative <- seConservative^2
  if(var_conservative - adjustmentFactor < 0){
    if(var_conservative != 0){
    base::warning("var_conservative is less than adjustmentFactor")
    }
    se_adjusted <- 0
  }else{
    se_adjusted <- base::sqrt( var_conservative - adjustmentFactor )
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
    t_less_than_g_index <- base::which(t_list < g)

    #If no periods less than g, there are no restrictions, return null
    if(base::length(t_less_than_g_index) == 0){
      return(NULL)
    }

    #Otherwise, Atilde is a matrix with rows= # periods before g and T columns.
    # It has 1s on the diagonal and 0s elsewhere. I.e. Atilde selects rows of Y corresponding with t<g
    Atilde0_g <- matrix(0,
                        nrow = base::length(t_less_than_g_index),
                        ncol = base::length(t_list))
    base::diag(Atilde0_g) <- 1

    return(Atilde0_g)
  }

  Atilde0_list <- purrr::map(.x = 1:(base::length(g_list)-1),
                             .f = createAtilde0_g)

  A0_list <- purrr::map(.x = 1:base::length(Atilde0_list),
                        .f = function(i){
                          #A0_g is a block matrix that stacks the Atilde_g's
                          #(for g < gmax) and multiplies by 0 for all g' \neq g
                          Atilde_times_indicator_list <- purrr::map(.x = 1:base::length(Atilde0_list),
                                                                    .f = ~(i == .x) * Atilde0_list[[.x]]  )
                          A0 <- purrr::reduce(Atilde_times_indicator_list,
                                              base::rbind)
                          return(A0)
                        }  )
  if(base::length(A0_list) == 1){
    A0_gmax <- -A0_list[[1]]}
  else{
    A0_gmax <- -base::Reduce(x = A0_list,
                             f = '+')
  }

  A0_list[[base::length(A0_list) + 1]] <- A0_gmax
  return(A0_list)
}


sum_of_lists <- function(.l){
  #Given a list of lists all of the same length, this returns a list where the jth entry is the sum of the jth entry of all the lists
  if(base::length(.l) == 1){return(.l[[1]])}
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
  # Given a list of lists all of the same length, where each inner list
  # is a vector, this returns a list where the jth element is a matrix
  # stacking all the jth elemtns
  results <- purrr::reduce(.x = .l,
                           .f = function(l1,
                                         l2){
                             purrr::map2(l1,
                                         l2,
                                         .f = ~base::rbind(.x , .y))
                           }
  )
  return(results)
}

# ----------------------------------------------------------------------
# Helpers for making the A_theta, A_0 lists
# ----------------------------------------------------------------------

create_Atheta_list_for_ATE_tg <- function(t,
                                          g,
                                          g_list,
                                          t_list,
                                          N_g_list,
                                          use_last_treated_only = FALSE,
                                          showWarnings = TRUE){

  #if(t < g & showWarnings){base::warning("t is less than g. ATE(t,g) is zero by assumption")}
  if(t >= base::max(g_list)){base::stop("t is greater than max(g)-1; ATE(t,g) is not identified.")}


  # Create A_thetas for ATT_{t,g}
  # -----------------------------

  treated_cohort_index <- base::which(g_list == g)
  if(!use_last_treated_only){
    #Create a list of which cohorts are eligible to be controls for each of the cohorts
    #This will be a null list if not eligible
    # To make placebos more comparable, we allow the control cohorts to be only those treated after max of g,t
    # This prevents cohort g from being include in its own control group when doing placebos
    control_cohort_indices <- base::which(g_list > base::pmax(g,t))
  }else{
    #If use_last_treated_only, compare only to the last treated cohort (i.e. max(G))
    control_cohort_indices <- base::which((g_list > t) & (g_list == base::max(g_list)))
  }
  N_control_cohorts <- base::sum(N_g_list[control_cohort_indices])

  #Weights control cohort
  A_theta_t <- numeric(base::length(g_list))
  A_theta_t[control_cohort_indices] <- - N_g_list[control_cohort_indices] / N_control_cohorts
  #Treated cohort
  A_theta_t[treated_cohort_index] <- A_theta_t[treated_cohort_index] + 1
  return(A_theta_t)
}


create_Atheta_list_for_event_study <- function(eventTime,
                                               g_list,
                                               t_list,
                                               N_g_DT,
                                               use_last_treated_only = FALSE){

  # Create A_thetas for an ``event-study'' coefficient at lag eventTime
  # This is the average treatment effects for units eventTime periods after first being treated
  # The average is taken over all cohorts g such that there is some untreated cohort at g+eventTime
  # Cohorts are weighted by the cohort size (N_g)

  maxG <- base::max(g_list)
  eligible_cohort_index <- base::which( ( pmax(g_list + eventTime, g_list) < maxG ) & ( pmax(g_list + eventTime, g_list) <= base::max(t_list) ) )
  N_g_list <- N_g_DT$N_g

  if(base::length(eligible_cohort_index) == 0){
    base::stop("There are no comparison cohorts for the given eventTime")
  }

  A_theta <- matrix(0, base::length(g_list), base::length(t_list))
  N_eligible <- base::sum(N_g_list[eligible_cohort_index])

  # length-G vector with weights for t, g, scaled by N_g/sum N_g
  create_Atheta_tg_helper <- function(gIndex) {
    g <- g_list[gIndex]
    (N_g_list[gIndex]/N_eligible) *
      create_Atheta_list_for_ATE_tg(g+eventTime, g, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # the results are placed in the time period columns corresponding to
  # the eligible cohorts offset by eventTime
  tindices <- base::unlist(base::sapply(g_list[eligible_cohort_index]+eventTime,
                                        function(t) base::which(t_list == t)))

  # we can mechanically run the helper function for all eligible cohorts,
  # but not  all eligible cohorts will have a corresponding time column
  # (this is particularly plausible when the event time offset is large)
  gindices <- base::unlist(base::sapply(eligible_cohort_index,
                                        function(g) if (base::any(t_list == g_list[g]+eventTime)) g else NULL))
  A_theta[,tindices] <- base::sapply(eligible_cohort_index, create_Atheta_tg_helper)[,gindices]

  return(A_theta)
}


create_Atheta_list_for_ATE_calendar_t <- function(t,
                                                  g_list,
                                                  t_list,
                                                  N_g_list,
                                                  use_last_treated_only = FALSE){
  treated_by_t_indices <- base::which(g_list <= t)
  N_total_treated <- base::sum(N_g_list[treated_by_t_indices])

  # length-G vector with weights for t, g, scaled by N_g/sum N_g
  create_Atheta_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_treated) *
          create_Atheta_list_for_ATE_tg(t, g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }

  # Here t is fixed, so all the G x 1 vectors correspond to the same
  # time column, and must therefore be aggregated
  A_theta_t <- base::Reduce(`+`, base::lapply(treated_by_t_indices, create_Atheta_tg_helper))

  # Returns single G x 1 vector, corresponding to t-th column
  return(A_theta_t)
}


create_Atheta_list_for_calendar_average_ATE <- function(g_list,
                                                        t_list,
                                                        N_g_DT,
                                                        use_last_treated_only = FALSE){

  t_eligible_index <- base::which((t_list >= base::min(g_list)) & (t_list < base::max(g_list)))
  T_eligible <- base::length(t_eligible_index)
  N_g_list <- N_g_DT$N_g
  A_theta <- base::matrix(0, base::length(g_list), base::length(t_list))

  # G x 1 vector containing the sum of all such vectors for
  # the t-th time period
  create_Atheta_tg_helper <- function(t) {
      create_Atheta_list_for_ATE_calendar_t(t, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # Place each vector into the column corresponding to each
  # time period, scaled by 1/T_eligible
  A_theta[,t_eligible_index] <-
    base::sapply(t_list[t_eligible_index], create_Atheta_tg_helper) / T_eligible

  return(A_theta)
}


create_Atheta_list_for_ATE_cohort_g <- function(g,
                                                g_list,
                                                t_list,
                                                N_g_list,
                                                use_last_treated_only = FALSE){
  treated_period_indices <- base::which(t_list >= g & t_list < base::max(g_list))
  T_treated <- base::length(treated_period_indices)
  A_theta_g <- base::matrix(0, base::length(g_list), base::length(t_list))

  # length-G vector with weights for t, g; will scale
  # later by 1/T_treated
  create_Atheta_tg_helper <- function(t) {
      create_Atheta_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # g is fixed; compute the vector for each treated period
  # and place in the corresponding columns
  A_theta_g[,treated_period_indices] <-
    base::sapply(t_list[treated_period_indices], create_Atheta_tg_helper)/T_treated

  # This is a matrix, which will be aggregated by the parent function
  return(A_theta_g)
}


create_Atheta_list_for_cohort_average_ATE <- function(g_list,
                                                      t_list,
                                                      N_g_DT,
                                                      use_last_treated_only = FALSE){
  g_eligible_index <-  base::which((g_list < base::max(g_list)) & (g_list <= base::max(t_list)))
  N_total_eligible <- N_g_DT[g_eligible_index, base::sum(N_g)]
  N_g_list <- N_g_DT$N_g

  # G x T matrix with all relevant time columns populated for the
  # g-th cohort, scaled by N_g / sum N_g
  create_Atheta_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_eligible) *
          create_Atheta_list_for_ATE_cohort_g(g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }

  # Aggregate G x T matrices across cohorts and return result
  A_theta <- base::Reduce(`+`, base::lapply(g_eligible_index, create_Atheta_tg_helper))
  return(A_theta)
}


create_Atheta_list_for_simple_average_ATE <- function(g_list,
                                                      t_list,
                                                      N_g_DT,
                                                      use_last_treated_only = FALSE){

  N_g_list <- N_g_DT$N_g
  #Create a df with all the (g,t) pairs for which ATE is identified;
  #Join in N_g for each of these pairs
  gt_df <- N_g_DT[data.table::CJ(g=g_list, t=t_list)[t >= g & t < base::max(g_list)]]

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- gt_df[,base::sum(N_g)]
  A_theta <- base::matrix(0, base::length(g_list), base::length(t_list))
  gt_df[,N_g := N_g / N_total]

  #sort by t
  data.table::setkey(gt_df, t, g)

  # generate the requisite G x 1 vector for each t, g pair, scaled by
  # N_g / sum N_g (note the variable N_g was replaced above)
  create_Atheta_tg_helper <- function(t, g, N_g) {
      N_g * create_Atheta_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # This creates a dataset for each t, g pair with G additional columns;
  # basically we've made each row of the dataset a 1 x G vector
  A_theta_gt <- gt_df[, base::as.list(create_Atheta_tg_helper(t, g, N_g)), by=c("t", "g")]
  data.table::setnames(A_theta_gt, c("t", "g", g_list))

  # If we sum the columns of this data set by t, we are summing up all
  # the vectors that correspond to that column
  A_theta_t <- A_theta_gt[,base::lapply(.SD, sum), by="t", .SDcols=base::as.character(g_list)]

  # finally, place the results of the sum into a G x T matrix;
  # the groups are the time periods, so we simply place the result
  # into the columns corresponding to those time periods
  indices <- base::unlist(base::sapply(A_theta_t$t, function(t) base::which(t_list == t)))
  A_theta[,indices] <- base::t(A_theta_t[,!"t"])

  return(A_theta)
}

# ----------------------------------------------------------------------
# Misc helpers
# ----------------------------------------------------------------------

calculate_full_vcv <- function(A_theta_list_list,
                               A_0_list_list,
                               Ybar_g_list,
                               S_g_list,
                               N_g_list,
                               g_list,
                               t_list,
                               beta,
                               resultsDF){

  # unless the user speficied a beta, compute beta for each A_theta,
  # A_0; note A_theta_list_list and A_0_list_list are named as such
  # because they are list of lists, since this code is meant to combine
  # multiple estimators (concretely, several eventstudy estimators).
  beta_list <-
  base::lapply(1:base::length(A_theta_list_list),
         function(i) {
           if(base::is.null(beta)) {
             compute_Betastar(Ybar_g_list,
                              A_theta_list_list[[i]],
                              A_0_list_list[[i]],
                              S_g_list,
                              N_g_list)
           }
           else {
             beta
           }
         })

  # Calculate A_theta_list - A_0_list %*% beta for each list
  # --------------------------------------------------------

  # The resulting list combined_A_list is a list of matrices of length
  # |G| so that the vector of thetas is \sum combined_A_list[g] %*% Ybar

  combine_A_lists <- function(A_theta_list, A_0_list, beta) {
    # Compute beta' %*% A_0_list
    A_0_beta_list <- left_product_lists(c = base::t(beta), .l = A_0_list)
    # Compute A_theta_list - A_0_list %*% beta
    combined_list <- sum_of_lists(list(A_theta_list, scalar_product_lists(-1, A_0_beta_list)) )
    return(combined_list)
  }
  combined_A_list <-
    stack_rows_of_lists(base::asplit(base::mapply(combine_A_lists, A_theta_list_list, A_0_list_list, beta_list), 2))

  # The newman vcv is \sum 1/N_g * combined_A * S * combined_A'
  vcv_neyman_terms_list <-
    purrr::pmap(.l = list(S_g = S_g_list, A = combined_A_list, N_g = base::as.list(N_g_list)),
                .f = function(S_g, A, N_g){ return((1/N_g) * A %*% S_g %*% base::t(A) ) } )
  vcv_neyman <- base::Reduce(f = '+', x = vcv_neyman_terms_list)

  # Next, we compute the refinement to the variance estimator
  # ---------------------------------------------------------

  # First, we find the earliest gmin across the different event-study coefs
  gMin <- base::min(base::mapply(function(x) computeGMin(x, g_list), A_theta_list_list))
  combine_betahat_lists <- function(A_theta_list, A_0_list, beta) {
    compute_se_Thetahat_beta(beta            = beta,
                             Ybar_g_list     = Ybar_g_list,
                             A_theta_list    = A_theta_list,
                             A_0_list        = A_0_list,
                             S_g_list        = S_g_list,
                             N_g_list        = N_g_list,
                             g_list          = g_list,
                             t_list          = t_list,
                             return_beta_sum = TRUE,
                             gMin            = gMin)
  }
  betahat_g_list <-
    base::asplit(base::mapply(combine_betahat_lists, A_theta_list_list, A_0_list_list, beta_list), 2)
  stacked_betahat_g_sum <- purrr::map(.x = betahat_g_list, .f = ~base::t(.x$betahat_g_sum)) %>%
                           base::Reduce(x = ., f = base::rbind)

  N <- base::sum(N_g_list)
  vcv_adjustment <- (1/N) * stacked_betahat_g_sum %*% betahat_g_list[[1]]$avg_MSM %*% base::t(stacked_betahat_g_sum)
  vcv <- vcv_neyman - vcv_adjustment

  #Set negative variances, if any, to 0 (arises from numerical precision issues sometimes)
  var_negative <- base::which( base::diag(vcv) < 0 )
  vcv[var_negative, var_negative] <- 0

  return(list(vcv = vcv, vcv_neyman = vcv_neyman))
}


processDF <- function(df, i, g, t, y){

  # Processes df inputted to staggered (or staggered_cs/sa)
  # -------------------------------------------------------
  #
  # - Checks the columns in the user-inputted values of i,g,t,y are actually in the data
  # - Rename the columns to "i", "g", "t", "y"
  # - Return a data.table for performance

  # Let's make sure we have columns with name i, t, y and g
  colnames_df <- base::colnames(df)
  if(!i %in% colnames_df){
    base::stop(base::paste0("There is no column ", i, " in the data. Thus, we are not able to find the unit identifier variable."))
  }
  if(!t %in% colnames_df){
    base::stop(base::paste0("There is no column ", t, " in the data. Thus, we are not able to find the time identifier variable."))
  }
  if(!g %in% colnames_df){
    base::stop(base::paste0("There is no column ", g, " in the data. Thus, we are not able to find the group identifier variable."))
  }
  if(!y %in% colnames_df){
    base::stop(base::paste0("There is no column ", y, " in the data. Thus, we are not able to find the outcome variable."))
  }

  # Sanity checks
  if(i %in% c("g", "t", "y" )){
    base::stop(base::paste0("Unit identifier cannot be labeled g, t, or y"))
  }

  if(t %in% c("i","y", "g" )){
    base::stop(base::paste0("Time identifier cannot be labeled i, g, or y"))
  }

  if(g %in% c("i", "t" ,"y" )){
    base::stop(base::paste0("Group identifier cannot be labeled i, t, or y"))
  }

  # Re-label i, t, g, y
  oldnames <- c(i, t, g, y)
  newnames <- c("i", "t", "g", "y")
  indices  <- base::which(!oldnames %in% newnames)
  oldnames <- oldnames[indices]
  newnames <- newnames[indices]
  if ( base::length(newnames) ) {
    dropnames <- base::intersect(newnames, base::names(df))
    if ( base::length(dropnames) ) {
      df <- df[,-base::match(dropnames, base::names(df))]
    }
    base::names(df)[base::match(oldnames, base::names(df))] <- newnames
  }

  # data.table for speed
  return(data.table::data.table(df)[,c("i", "t", "g", "y")])
}

# ----------------------------------------------------------------------
# Main functions
# ----------------------------------------------------------------------

#' @useDynLib staggered
#' @importFrom data.table ":="
#' @importFrom magrittr "%>%"
#' @import Rcpp
#' @import RcppEigen
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
#' @param use_last_treated_only If true, then A_0_list and A_theta_list are created to only make comparisons with the last treated cohorts (as suggested by Sun and Abraham), rather than using not-yet-treated units as comparisons. If set to TRUE (and use_DiD_A0 = TRUE), then beta=1 corresponds with the Sun and Abraham estimator.
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.
#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman).
#' @references
#' \cite{Roth, Jonatahan, and Sant'Anna, Pedro H. C. (2021),
#'   'Efficient Estimation for Staggered Rollout Designs', arXiv: 2102.01291, \url{https://arxiv.org/abs/2102.01291}.}
#' @examples
#' \dontshow{
#'   # restrict threads for CRAN compliance
#'   dt_threads <- data.table::getDTthreads()
#'   data.table::setDTthreads(1)
#' }
#'
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
#' head(eventPlotResults)
#' \dontshow{
#'   # restore thread setting
#'   data.table::setDTthreads(dt_threads)
#' }
#'
#' @export
staggered <- function(df,
                      i = "i",
                      t = "t",
                      g = "g",
                      y = "y",
                      estimand     = NULL,
                      A_theta_list = NULL,
                      A_0_list     = NULL,
                      eventTime    = 0,
                      beta         = NULL,
                      use_DiD_A0   = ifelse(is.null(A_0_list), TRUE, FALSE),
                      return_full_vcv         = FALSE,
                      use_last_treated_only   = FALSE,
                      compute_fisher          = FALSE,
                      num_fisher_permutations = 500,
                      skip_data_check         = FALSE){

  stats::runif(1) # make sure .Random.seed exists
  rseed.cached <- .Random.seed
  base::on.exit({.Random.seed <<- rseed.cached})

  # If estimand is provided, force to be lower-case (allowing for
  # non-case sensitive inputs)
  if(!base::is.null(estimand)){
    estimand <- base::tolower(estimand)
  }

  if(base::length(eventTime) > 1) {
    if(estimand != "eventstudy"){
      base::stop("You provided a vector fpr eventTime but estimand is not set to 'eventstudy'. Did you mean to set estimand = 'eventstudy'?")
    }
  }

  #Save the user-inputted beta (used in FRT and event study calls)
  user_input_beta <- beta

  #-----------------------
  #Process the inputted df
  #-----------------------

  #Check inputted columns and renaming to i,t,g,y, and balancing on
  #(i,t); skip if skip_data_check = TRUE
  if(!skip_data_check){
    df <- processDF(df,
                    i=i,
                    g=g,
                    t=t,
                    y=y)
    if ( compute_fisher ) {
      df$cached_order <- 1:base::NROW(df)
    }
    #Balance the panel (and throw a warning if original panel is unbalanced)
    df <- balance_df(df = df)
  } else {
    if ( compute_fisher ) {
      df$cached_order <- 1:base::NROW(df)
    }
  }

  # --------------------------
  # Flag for singleton cohorts
  # --------------------------

  # Compute number of units per cohort
  cohort_size <- df[,.(N=data.table::uniqueN(i)),by="g"]
  flag_singleton <- cohort_size[N == 1,g]
  l_flag1 <- base::length(flag_singleton)
  # Drop cohorts which are singleton
  if ( l_flag1 ) {
    gpaste <- base::paste(flag_singleton, collapse=", ")
    if(l_flag1==1){
      base::warning(base::paste0("The treatment cohort g = ", gpaste, " has a single cross-sectional unit. We drop this cohort."))
    } else {
      base::warning(base::paste0("The treatment cohorts g = ", gpaste, " have a single cross-sectional unit only. We drop these cohorts."))
    }
    df <- df[!(g %in% flag_singleton)]
  }

  # -----------------------
  # Set up helper variables
  # -----------------------

  g_level_summaries <- compute_g_level_summaries(df, is_balanced = TRUE)
  Ybar_g_list <- g_level_summaries$Ybar_g_List
  S_g_list    <- g_level_summaries$S_g_List
  N_g_DT      <- g_level_summaries$N_g_DT
  g_list      <- g_level_summaries$g_list
  t_list      <- g_level_summaries$t_list
  N_g_list    <- N_g_DT$N_g

  #If estimand is provided, calculate the appropriate A_theta_list
  if(!base::is.null(estimand)){
    if(estimand == "simple"){
      A_theta <- create_Atheta_list_for_simple_average_ATE(g_list = g_list,
                                                           t_list = t_list,
                                                           N_g_DT = N_g_DT,
                                                           use_last_treated_only = use_last_treated_only)
      A_theta_list <- base::lapply(base::asplit(A_theta, 1), base::rbind)
    }else if(estimand == "cohort"){
      A_theta <- create_Atheta_list_for_cohort_average_ATE(g_list = g_list,
                                                           t_list = t_list,
                                                           N_g_DT = N_g_DT,
                                                           use_last_treated_only = use_last_treated_only)
      A_theta_list <- base::lapply(base::asplit(A_theta, 1), base::rbind)
    }else if(estimand == "calendar"){
      A_theta <- create_Atheta_list_for_calendar_average_ATE(g_list = g_list,
                                                             t_list = t_list,
                                                             N_g_DT = N_g_DT,
                                                             use_last_treated_only = use_last_treated_only)
      A_theta_list <- base::lapply(base::asplit(A_theta, 1), base::rbind)
    }else if(estimand == "eventstudy"){
      create_Atheta_es_helper <- function(eventTime) {
          A_theta <-
          create_Atheta_list_for_event_study(eventTime = eventTime,
                                             g_list = g_list,
                                             t_list = t_list,
                                             N_g_DT = N_g_DT,
                                             use_last_treated_only = use_last_treated_only)
          base::lapply(base::asplit(A_theta, 1), base::rbind)
      }
      A_theta_list <- if ( base::length(eventTime) > 1 ) {
        base::lapply(eventTime, create_Atheta_es_helper)
      } else {
        create_Atheta_es_helper(eventTime)
      }
    }
  }

  #If no valid estimand is provided and no A_theta_list, throw and error
  if(base::is.null(A_theta_list)){
    base::stop("Estimand must be one of simple, cohort, calendar, or eventstudy; or custom A_theta_list must be provided")
  }

  #Create A_0_list if a custom A_0_list is not provided
  if(base::is.null(A_0_list) & (use_DiD_A0==FALSE)){
    A_0_list <- create_A0_list(g_list = g_list, t_list = t_list)
  }

  #If use_DiD_A0, use only the A0's associated with the DiD estimand
  if(use_DiD_A0){

    if(base::is.null(estimand)){
      base::stop("If use_DiD_A0 = TRUE, you must provide an estimand.")
    }

    if(estimand == "simple"){
      A_0 <- create_A0_list_for_simple_average_ATE(g_list = g_list,
                                                   t_list = t_list,
                                                   N_g_DT = N_g_DT,
                                                   use_last_treated_only = use_last_treated_only)
      A_0_list <- base::lapply(base::asplit(A_0, 1), base::rbind)
    }else if(estimand == "cohort"){
      A_0 <- create_A0_list_for_cohort_average_ATE(g_list = g_list,
                                                   t_list = t_list,
                                                   N_g_DT = N_g_DT,
                                                   use_last_treated_only = use_last_treated_only)
      A_0_list <- base::lapply(base::asplit(A_0, 1), base::rbind)
    }else if(estimand == "calendar"){
      A_0 <- create_A0_list_for_calendar_average_ATE(g_list = g_list,
                                                     t_list = t_list,
                                                     N_g_DT = N_g_DT,
                                                     use_last_treated_only = use_last_treated_only)
      A_0_list <- base::lapply(base::asplit(A_0, 1), base::rbind)
    }else if(estimand == "eventstudy"){
      create_A0_es_helper <- function(eventTime) {
          A_0 <-
          create_A0_list_for_event_study(eventTime = eventTime,
                                         g_list = g_list,
                                         t_list = t_list,
                                         N_g_DT = N_g_DT,
                                         use_last_treated_only = use_last_treated_only)
          base::lapply(base::asplit(A_0, 1), base::rbind)
      }
      A_0_list <- if ( base::length(eventTime) > 1 ) {
        base::lapply(eventTime, create_A0_es_helper)
      } else {
        create_A0_es_helper(eventTime)
      }
    }
  }

  # ----------
  # FRT helper
  # ----------

  # The point of this helper function is to loop over the number of
  # fisher repetitions; it calls staggered recursively, so it's slow
  compute_fisher_fun <- function(A_theta_list, A_0_list, eventTime, i_g_table, seed) {
    # This function takes a data.frame with columns i and g, and
    # permutes the values of g assigned to i
    #
    # The input i_g_table has the unique combinations of (i,g) in df,
    # and is calculated outside for speed improvements

    # Draw a random permutation of the elements of first_period_df
    # Replace first_period_df$g with a permuted version based on randIndex
    base::set.seed(seed)
    i_g_table[, g := g_cached[base::sample(.N)]]

    # Merge the new treatment assignments back with the original
    # (NB: df is modified in place)
    df[i_g_table[,!"g_cached"], g := (i.g)]

    staggered::staggered(df                 = df,
                         estimand           = NULL,
                         beta               = user_input_beta,
                         A_theta_list       = A_theta_list,
                         A_0_list           = A_0_list,
                         eventTime          = eventTime,
                         return_full_vcv    = FALSE,
                         compute_fisher     = FALSE,
                         skip_data_check    = TRUE)
  }

  # --------------------
  # Final results helper
  # --------------------

  # The main reason to use this helper function is to loop over
  # the number of event times; note that while it takes eventTime
  # as argument for consistency, strictly speaking it shouldn't be
  # necessary, since everything has been pre-computed
  resultsFun <- function(A_theta_list, A_0_list, eventTime) {
    Xvar_list <- purrr::pmap(.l = list(A_0_list, S_g_list, base::as.list(N_g_list)) ,
                             .f = function(A0,S,N){
                               return(1/N * eigenMapMatMult(eigenMapMatMult(A0,S), base::t(A0)))
                             })

    if(base::is.null(user_input_beta)){
      beta <- compute_Betastar(Ybar_g_list,
                               A_theta_list,
                               A_0_list,
                               S_g_list,
                               N_g_list,
                               Xvar_list)
    }

    #If beta =0, convert beta to the appropriate length
    if ( base::length(beta) == 1 && beta == 0){
      beta <- base::matrix(0, base::length(A_0_list))
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
                                          return_beta_sum = TRUE)

    se <- seResults$se

    resultsDF <- base::data.frame(estimate  = thetahat,
                                  se        = se,
                                  se_neyman = se_neyman)

    if ( compute_fisher ) {

      # Find unique pairs of (i,g). This will be used for computing the permutations
      i_g_table <- df[t == base::min(t), c("i", "g", "cached_order")]
      i_g_table[, c("g_cached", "cached_order") := list(g[order(cached_order)], NULL)]

      FRTResults <-
        purrr::map(.x = 1:num_fisher_permutations,
                   .f = purrr::possibly(
                          .f = ~compute_fisher_fun(A_theta_list,
                                                   A_0_list,
                                                   eventTime,
                                                   i_g_table,
                                                   .x),
                          otherwise = NULL)
        )
      FRTResults <- base::as.data.frame(base::do.call(base::rbind, FRTResults))

      # Now, we compute the FRT for each seed, permuting treatment for each one
      #
      # We catch any errors in the FRT simulations, and throw a warning if
      # at least one has an error (using the remaining draws to calculate frt)
      successful_frt_draws <- base::NROW(FRTResults)
      if (successful_frt_draws < num_fisher_permutations) {
        base::warning("There was an error in at least one of the FRT simulations. Removing the problematic draws.")
      }

      resultsDF$fisher_pval <- base::mean( base::abs(resultsDF$estimate/resultsDF$se) < base::abs(FRTResults$estimate/FRTResults$se) )
      resultsDF$fisher_pval_se_neyman <- base::mean( base::abs(resultsDF$estimate/resultsDF$se_neyman) < base::abs(FRTResults$estimate/FRTResults$se_neyman) )
      resultsDF$num_fisher_permutations <- successful_frt_draws
    }

    return(resultsDF)
  }

  #If eventTime is a vector, call staggered for each event-time and
  #combine the results Add the variable eventTime to the data frame
  if ( base::length(eventTime) > 1) {
    resultsDF <- base::do.call(base::rbind,
                               base::mapply(resultsFun,
                                            A_theta_list,
                                            A_0_list,
                                            base::as.list(eventTime),
                                            SIMPLIFY=FALSE))
    #Add in eventTimes
    resultsDF$eventTime <- eventTime
  } else {
    resultsDF <- resultsFun(A_theta_list, A_0_list, eventTime)
  }

  # either add the full vcv or return DF as is
  if( !return_full_vcv ) {
    return(resultsDF)
  } else {
    if(base::length(eventTime) > 1) {
      vcvs <- calculate_full_vcv(A_theta_list, A_0_list,
                                 beta        = user_input_beta,
                                 Ybar_g_list = Ybar_g_list,
                                 S_g_list    = S_g_list,
                                 N_g_list    = N_g_list,
                                 g_list      = g_list,
                                 t_list      = t_list,
                                 resultsDF   = resultsDF)

      resultsList <- list(resultsDF = resultsDF,
                          vcv = vcvs$vcv,
                          vcv_neyman = vcvs$vcv_neyman)
    }
    else {
      resultsList <- list(resultsDF  = resultsDF,
                          vcv        = base::as.matrix(resultsDF$se^2),
                          vcv_neyman = base::as.matrix(resultsDF$se_neyman^2))
    }
    return(resultsList)
  }
}


#' @title Calculate the Callaway & Sant'Anna (2021) estimator for staggered rollouts
#' @description This functions calculates the Callaway & Sant'Anna (2021) estimator for staggered rollout designs using not-yet-treated units (including never-treated, if available) as controls.
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
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.

#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman).
#' @references
#'   \cite{Callaway, Brantly, and Sant'Anna, Pedro H. C. (2021),
#'   'Difference-in-Differences with Multiple Time Periods', Journal of Econometrics,
#'   \doi{10.1016/j.jeconom.2020.12.001}.}
#' @examples
#' \dontshow{
#'   # restrict threads for CRAN compliance
#'   dt_threads <- data.table::getDTthreads()
#'   data.table::setDTthreads(1)
#' }
#'
#' # Load some libraries
#' set.seed(1234)
#' # load the officer data and subset it
#' df <- pj_officer_level_balanced
#' group_random <- sample(unique(df$assigned), 3)
#' df <- df[df$assigned %in% group_random,]
#' # We modify the data so that the time dimension is named t,
#' # the period of treatment is named g,
#' # the outcome is named y,
#' # and the individual identifiers are named i
#' # (this allow us to use default arguments on \code{staggered_cs}).
#' oldnames <- c("period", "complaints", "first_trained", "uid")
#' names(df)[match(oldnames, names(df))] <- c("t", "y", "g", "i")
#' # Calculate Callaway and Sant'Anna estimator for the simple weighted average
#' staggered_cs(df = df, estimand = "simple")
#' # Calculate Callaway and Sant'Anna estimator for the cohort weighted average
#' staggered_cs(df = df, estimand = "cohort")
#' # Calculate Callaway and Sant'Anna estimator for the calendar weighted average
#' staggered_cs(df = df, estimand = "calendar")
#' # Calculate Callaway and Sant'Anna event-study coefficients for the first 24 months
#' # (month 0 is instantaneous effect)
#' eventPlotResults <- staggered_cs(df = df, estimand = "eventstudy", eventTime = 0:23)
#' head(eventPlotResults)
#' \dontshow{
#'   # restore thread setting
#'   data.table::setDTthreads(dt_threads)
#' }
#'
#' @export
staggered_cs <- function(df,
                         i = "i",
                         t = "t",
                         g = "g",
                         y = "y",
                         estimand     = NULL,
                         A_theta_list = NULL,
                         A_0_list     = NULL,
                         eventTime    = 0,
                         return_full_vcv         = FALSE,
                         compute_fisher          = FALSE,
                         num_fisher_permutations = 500,
                         skip_data_check         = FALSE){

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
  TreatedBeforeMinT <- df[g <= base::min(t), .N]
  if( TreatedBeforeMinT ) {
    df <- df[g > base::min(t)]
    warning("Dropping units who were treated in the first period or earlier, since CS estimator is not defined (and ATT(t,g) not identified under parallel trends).")
  }

  results <- staggered(df                      = df,
                       estimand                = estimand,
                       A_theta_list            = A_theta_list,
                       eventTime               = eventTime,
                       beta                    = 1,
                       use_DiD_A0              = TRUE,
                       use_last_treated_only   = FALSE,
                       return_full_vcv         = return_full_vcv,
                       compute_fisher          = compute_fisher,
                       num_fisher_permutations = num_fisher_permutations,
                       skip_data_check         = skip_data_check)

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
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.
#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman).
#' @references
#'   \cite{Sun, Liyang, and Abraham, Sarah (2020),
#'   'Estimating dynamic treatment effects in event studies with heterogeneous treatment effects', Forthcoming at the Journal of Econometrics,
#'   \doi{10.1016/j.jeconom.2020.09.006}.}
#' @examples
#' \dontshow{
#'   # restrict threads for CRAN compliance
#'   dt_threads <- data.table::getDTthreads()
#'   data.table::setDTthreads(1)
#' }
#'
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
#' oldnames <- c("period", "complaints", "first_trained", "uid")
#' names(df)[match(oldnames, names(df))] <- c("t", "y", "g", "i")
#' # Calculate Sun and Abraham estimator for the simple weighted average
#' staggered_sa(df = df, estimand = "simple")
#' # Calculate Sun and Abraham estimator for the cohort weighted average
#' staggered_sa(df = df, estimand = "cohort")
#' # Calculate Sun and Abraham estimator for the calendar weighted average
#' staggered_sa(df = df, estimand = "calendar")
#' # Calculate Sun and Abraham event-study coefficients for the first 24 months
#' # (month 0 is instantaneous effect)
#' eventPlotResults <- staggered_sa(df = df, estimand = "eventstudy", eventTime = 0:23)
#' head(eventPlotResults)
#' \dontshow{
#'   # restore thread setting
#'   data.table::setDTthreads(dt_threads)
#' }
#'
#' @export
staggered_sa <- function(df,
                         i = "i",
                         t = "t",
                         g = "g",
                         y = "y",
                         estimand     = NULL,
                         A_theta_list = NULL,
                         A_0_list     = NULL,
                         eventTime    = 0,
                         return_full_vcv         = FALSE,
                         compute_fisher          = FALSE,
                         num_fisher_permutations = 500,
                         skip_data_check         = FALSE){


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
  TreatedBeforeMinT <- df[g <= base::min(t), .N]
  if( TreatedBeforeMinT ) {
    df <- df[g > base::min(t)]
    warning("Dropping units who were treated in the first period or earlier, since SA estimator is not defined (and ATT(t,g) not identified under parallel trends).")
  }

  results <- staggered(df                      = df,
                       estimand                = estimand,
                       A_theta_list            = A_theta_list,
                       eventTime               = eventTime,
                       beta                    = 1,
                       use_DiD_A0              = TRUE,
                       use_last_treated_only   = TRUE,
                       return_full_vcv         = return_full_vcv,
                       compute_fisher          = compute_fisher,
                       num_fisher_permutations = num_fisher_permutations,
                       skip_data_check         = skip_data_check)

  return(results)
}
