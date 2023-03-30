
#' @title Calculate balance checks in staggered rollout designs
#' @description This functions calculates Wald-tests for balance in staggered rollout designs proposed by Roth and Sant'Anna.
#' @param df A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)
#' @param i The name of column containing the individual (cross-sectional unit) identifier. Default is "i".
#' @param t The name of the column containing the time periods. Default is "t".
#' @param g The name of the column containing the first period when a particular observation is treated, with Inf denoting never treated. Default is "g".
#' @param y The name of the column containing the outcome variable. Default is "y".
#' @param estimand The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. "eventstudy" returns the average effect at the ''event-time'' given in the parameter EventTime.  The parameter can be left blank if a custom parameter is provided in A_theta_list. The argument is not case-sensitive.
#' @param A_0_list This parameter allow for specifying the matrices used to construct the Xhat vector of pre-treatment differences. If left NULL, the default is to use the scalar set of controls used in Callaway and Sant'Anna. If use_DiD_A0 = FALSE, then it uses the full vector possible comparisons of (g,g') in periods t<g,g'.
#' @param eventTime If using estimand = "eventstudy", specify what eventTime you want the event-study parameter for. The default is 0, the period in which treatment occurs. If a vector is provided, estimates are returned for all the event-times in the vector.
#' @param use_DiD_A0 If this parameter is true, then Xhat corresponds with the scalar used by Callaway and Sant'Anna, so the Callaway and Sant'Anna estimator corresponds with beta=1. If it is false, the Xhat is a vector with all possible comparisons of pairs of cohorts before either is treated. The latter option should only be used when the number of possible comparisons is small relative to sample size.
#' @param use_last_treated_only If true, then A_0_list and A_theta_list are created to only make comparisons with the last treated cohorts (as suggested by Sun and Abraham), rather than using not-yet-treated units as comparisons. If set to TRUE (and use_DiD_A0 = TRUE), then beta=1 corresponds with the Sun and Abraham estimator.
#' @param compute_fisher If true, computes a Fisher Randomization Test using the studentized estimator.
#' @param num_fisher_permutations The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.
#' @param return_full_vcv If this is true, then the function returns a list containing the full variance-covariance matrix for all Xhats.
#' @param skip_data_check If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.
#' @param seed Set seed for permutations
#' @return resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman). (If return_matrix_list = TRUE, it likewise returns a list containing lists of matrices used in the vcv calculation.)
#' @references
#' \cite{Roth, Jonatahan, and Sant'Anna, Pedro H. C. (2021),
#'   'Efficient Estimation for Staggered Rollout Designs', arXiv: 2102.01291, \url{https://arxiv.org/abs/2102.01291}.}
#' @examples
#' set.seed(1234)
#' # load the officer data and subset it
#' df <- pj_officer_level_balanced
#' group_random <- sample(unique(df$assigned), 3)
#' df <- df[df$assigned %in% group_random,]
#' # Calculate balance checks for simple aggregation
#' balance_checks(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "simple")
#' # Calculate balance checks for the cohort weighted average
#' balance_checks(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "cohort")
#' # Calculate balance checks for the calendar weighted average
#' balance_checks(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "calendar")
#' # Calculate balance checks for event-study aggregation for the first 24 months
#' # (month 0 is instantaneous effect)
#' event_bal_checks <- balance_checks(df = df,
#'   i = "uid",
#'   t = "period",
#'   g = "first_trained",
#'   y = "complaints",
#'   estimand = "eventstudy",
#'   eventTime = 0:23)
#' head(event_bal_checks)
#'
#' @export
balance_checks <- function(df,
                           i = "i",
                           t = "t",
                           g = "g",
                           y = "y",
                           estimand = NULL,
                           A_0_list = NULL,
                           eventTime = 0,
                           use_DiD_A0 = NULL,
                           use_last_treated_only = FALSE,
                           compute_fisher = FALSE,
                           num_fisher_permutations = 500,
                           return_full_vcv = FALSE,
                           skip_data_check = FALSE,
                           seed = NULL){

  if(base::is.null(use_DiD_A0)) use_DiD_A0 <- base::ifelse(base::is.null(A_0_list), TRUE, FALSE)

  # If estimand is provided, force to be lower-case (allowing for
  # non-case sensitive inputs)
  if(!base::is.null(estimand)){
    estimand <- base::tolower(estimand)
  }

  # If eventTime is a vector, call staggered for each event-time and
  # combine the results Add the variable eventTime to the data frame
  if(base::length(eventTime) > 1){
    if(estimand != "eventstudy" & estimand != "all" ){
      base::stop("You provided a vector for eventTime but estimand is not set to 'eventstudy' or 'all'. Did you mean to set estimand = 'eventstudy' or estimand = 'all'?")
    }
  }

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
    #Balance the panel (and throw a warning if original panel is unbalanced)
    df <- balance_df(df = df)
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

  #Create A_0_list if a custom A_0_list is not provided
  if(base::is.null(A_0_list) & (use_DiD_A0==FALSE)){
    A_0_list <- create_A0_list(g_list = g_list, t_list = t_list)
  }

  #If use_DiD_A0, use only the A0's associated with the DiD estimand
  #(NB: Any number of estimands can be requested in this case, so we
  #append all the requisite matrices.)
  if(use_DiD_A0){

    if(base::is.null(estimand)){
      base::stop("If use_DiD_A0 = TRUE, you must provide an estimand.")
    }

    estim <- c()
    A_0_list <- list()
    if(base::length(base::intersect(estimand, c("simple", "all")))){
      A_0 <- create_A0_list_for_simple_average_ATE(g_list = g_list,
                                                   t_list = t_list,
                                                   N_g_DT = N_g_DT,
                                                   use_last_treated_only = use_last_treated_only)
      A_0_list <- c(A_0_list, list(base::lapply(base::asplit(A_0, 1), base::rbind)))
      estim <- c(estim, "all_simple")
    }
    if(base::length(base::intersect(estimand, c("cohort", "all")))){
      A_0 <- create_A0_list_for_cohort_average_ATE(g_list = g_list,
                                                   t_list = t_list,
                                                   N_g_DT = N_g_DT,
                                                   use_last_treated_only = use_last_treated_only)
      A_0_list <- c(A_0_list, list(base::lapply(base::asplit(A_0, 1), base::rbind)))
      estim <- c(estim, "all_cohort")
    }
    if(base::length(base::intersect(estimand, c("calendar", "all")))){
      A_0 <- create_A0_list_for_calendar_average_ATE(g_list = g_list,
                                                     t_list = t_list,
                                                     N_g_DT = N_g_DT,
                                                     use_last_treated_only = use_last_treated_only)
      A_0_list <- c(A_0_list, list(base::lapply(base::asplit(A_0, 1), base::rbind)))
      estim <- c(estim, "all_calendar")
    }
    if(base::length(base::intersect(estimand, c("eventstudy", "all")))){
      create_A0_es_helper <- function(eventTime) {
          A_0 <-
          create_A0_list_for_event_study(eventTime = eventTime,
                                         g_list = g_list,
                                         t_list = t_list,
                                         N_g_DT = N_g_DT,
                                         use_last_treated_only = use_last_treated_only)
          base::lapply(base::asplit(A_0, 1), base::rbind)
      }
      A_0_list <- c(A_0_list, base::lapply(eventTime, create_A0_es_helper))
      estim <- c(estim, base::paste0("all_ES",eventTime))
    }
  }

  # ------------------
  # Main balance check
  # ------------------

  A_0_list <- stack_rows_of_lists(A_0_list)
  balance_checks_Xhat <- compute_balance_test(Ybar_g_list = Ybar_g_list,
                                              A_0_list = A_0_list,
                                              S_g_list = S_g_list,
                                              N_g_list = N_g_list)
  se_Xhat <- base::sqrt(base::diag(base::as.matrix(balance_checks_Xhat$Xvar)))
  resultsDF <- data.frame(Xhat = balance_checks_Xhat$Xhat,
                          se_Xhat = se_Xhat,
                          t_test = base::abs(balance_checks_Xhat$Xhat/se_Xhat),
                          pvalue_t = 2*stats::pnorm(-base::abs(balance_checks_Xhat$Xhat/se_Xhat)),
                          Wald_test_Xhat = balance_checks_Xhat$Wald_test_Xhat,
                          pvalue_Wald = stats::pchisq(balance_checks_Xhat$Wald_test_Xhat,
                                                      df = base::length( balance_checks_Xhat$Xhat),
                                                      lower.tail = FALSE),
                          N = balance_checks_Xhat$N,
                          fisher_pval = NA,
                          fisher_supt_pval = NA,
                          num_fisher_permutations = NA,
                          estimand = estim)

  # --------------------
  # Do FRT, if specified
  # --------------------

  permuteTreatment2 <- function(df, i_g_table, seed = NULL) {
    # This function takes a data.frame with columns i and g, and
    # permutes the values of g assigned to i
    #
    # The input i_g_table has the unique combinations of (i,g) in df,
    # and is calculated outside for speed improvements

    # Draw a random permutation of the elements of first_period_df
    if(!base::is.null(seed)) base::set.seed(seed)
    i_g_table[, g := g[base::sample(.N)]]

    # Merge the new treatment assignments back with the original
    # (NB: df is modified in place)
    df[i_g_table, g := list(i.g)]
  }

  permutation_t_test <- function(df, i_g_table, A_0_list, seed = NULL){
    permuteTreatment2(df, i_g_table, seed)
    g_level_summaries <- compute_g_level_summaries(df, is_balanced = TRUE)
    Ybar_g_list <- g_level_summaries$Ybar_g_List
    S_g_list    <- g_level_summaries$S_g_List
    N_g_DT      <- g_level_summaries$N_g_DT
    g_list      <- g_level_summaries$g_list
    t_list      <- g_level_summaries$t_list
    N_g_list    <- N_g_DT$N_g

    balance_checks_Xhat <- compute_balance_test(Ybar_g_list = Ybar_g_list,
                                                A_0_list = A_0_list,
                                                S_g_list = S_g_list,
                                                N_g_list = N_g_list)
    se_Xhat <- base::sqrt(base::diag(base::as.matrix(balance_checks_Xhat$Xvar)))
    t_test <- base::abs(balance_checks_Xhat$Xhat/se_Xhat)
    matrix(t_test, nrow = 1)
  }

  FRTResults_bal = NULL
  if(compute_fisher == TRUE){
    # Find unique pairs of (i,g). This will be used for computing the permutations
    i_g_table <- df[t == base::min(t), c("i","g")]

    # Now, we compute the FRT for each seed, permuting treatment for each
    #
    # We catch any errors in the FRT simulations, and throw a warning if at
    # least one has an error (using the remaining draws to calculate frt)
    if(!base::is.null(seed)) base::set.seed(seed)
    seed_frt <-  1:num_fisher_permutations * base::floor(stats::rexp(1, rate = 1/500))

    FRTResults_bal <-
      purrr::map(.x = seed_frt,
                 .f = purrr::possibly(.f = ~permutation_t_test(df, i_g_table, A_0_list, .x),
                                      otherwise = NULL))
    FRTResults_bal <- base::do.call(base::rbind, FRTResults_bal)
    successful_frt_draws <- base::nrow(FRTResults_bal)

    if(successful_frt_draws < num_fisher_permutations){
      base::warning("There was an error in at least one of the FRT simulations. Removing the problematic draws.")
    }

    # Comput sup-t tests
    FRT_t <- base::as.numeric(base::apply(FRTResults_bal, 1, base::max))
    max_t <- base::max(resultsDF$t_test)
    resultsDF$fisher_supt_pval <- base::mean( max_t < FRT_t )

    # Need to compute p-value per row
    fisher_pval <- base::apply(FRTResults_bal, 1, '>', resultsDF$t_test)
    if(base::is.matrix(fisher_pval)) {
      fisher_pval <- base::rowMeans(fisher_pval)
    } else {
      fisher_pval <- base::mean(fisher_pval)
    }
    resultsDF$fisher_pval <- fisher_pval

    resultsDF$num_fisher_permutations <- successful_frt_draws
  }

  resultsDF <- base::as.data.frame(resultsDF)
  if(return_full_vcv){
    Xvar1 = balance_checks_Xhat$Xvar
  }
  else {
    Xvar1 = NULL
  }
  list_results <- list(resultsDF = resultsDF, Xvar = Xvar1, FRTResults = FRTResults_bal)

  return(list_results)
}
