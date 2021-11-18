
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
#' event_bal_checks %>% head()
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

  if(is.null(use_DiD_A0)) use_DiD_A0 <- ifelse(is.null(A_0_list),    TRUE,         FALSE)
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
  df_processed <- df


  #  Compute number of units per cohort
  cohort_size <- base::table(df_processed$g)/base::length(base::table(df_processed$t))
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

    df_processed <- df_processed[(df_processed$g %in% flag_singleton) == FALSE,]
  }


  #  If estimand is provided, force to be lower-case (allowing for non-case sensitive inputs)
  if(!is.null(estimand)){
    estimand <- tolower(estimand)
  }

  # #If eventTime is a vector, call staggered for each event-time and combine the results
  # #Add the variable eventTime to the data frame
  if(length(eventTime) > 1){

    if(estimand != "eventstudy" & estimand != "all" ){
      stop("You provided a vector for eventTime but estimand is not set to 'eventstudy' or 'all'. Did you mean to set estimand = 'eventstudy' or estimand = 'all'?")
    }
  }
  #
  #   eventPlotResultsList <-
  #     purrr::map(.x = eventTime,
  #                .f = ~balance_checks(df = df,
  #                                     estimand = estimand,
  #                                     A_0_list = A_0_list,
  #                                     eventTime = .x,
  #                                     use_DiD_A0 = use_DiD_A0,
  #                                     use_last_treated_only = use_last_treated_only,
  #                                     skip_data_check = TRUE))
  #
  #   resultsDF <- purrr::reduce(.x = purrr::map(.x = eventPlotResultsList, .f = ~ .x$resultsDF),
  #                              .f = dplyr::bind_rows)
  #
  #   #Add in eventTimes
  #   resultsDF$eventTime <- eventTime
  #
  #   return(resultsDF)
  # }

  g_level_summaries <- compute_g_level_summaries(df_processed, is_balanced = TRUE)
  Ybar_g_list <- g_level_summaries$Ybar_g_List
  S_g_list <- g_level_summaries$S_g_List
  N_g_list <- g_level_summaries$N_g_List
  g_list <- g_level_summaries$g_list
  t_list <- g_level_summaries$t_list




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

      estim <- "simple"
    }else if(estimand == "cohort"){
      A_0_list <- create_A0_list_for_cohort_average_ATE(g_list = g_list,
                                                        t_list = t_list,
                                                        N_g_list = N_g_list,
                                                        use_last_treated_only = use_last_treated_only)
      estim <- "cohort"
    }else if(estimand == "calendar"){
      A_0_list <- create_A0_list_for_calendar_average_ATE(g_list = g_list,
                                                          t_list = t_list,
                                                          N_g_list = N_g_list,
                                                          use_last_treated_only = use_last_treated_only)

      estim <- "calendar"
    }else if(estimand == "eventstudy"){
      # Event Study one
      if(length(eventTime)>1){
        n_eventTime = length(eventTime)
        aux_list = NULL
        for(j in 1:n_eventTime){
          aux_list_new <- create_A0_list_for_event_study(eventTime = eventTime[j],
                                                         g_list = g_list,
                                                         t_list = t_list,
                                                         N_g_list = N_g_list,
                                                         use_last_treated_only = use_last_treated_only)
          aux_list <- rbind(aux_list, aux_list_new)
        }

        A_0_list_event_study <- apply(aux_list, 2, rbind)

        A_0_list <-  purrr::pmap(.l = list(A_0_list_event_study),
                                 .f = function(A0){ return(do.call(rbind, A0)) }
        )

        estim <- paste0("ES",eventTime)


      } else {
        A_0_list <- create_A0_list_for_event_study(eventTime = eventTime,
                                                   g_list = g_list,
                                                   t_list = t_list,
                                                   N_g_list = N_g_list,
                                                   use_last_treated_only = use_last_treated_only)
        estim <- paste0("ES",eventTime)
      }

      #
      # A_0_list <- create_A0_list_for_event_study(eventTime = eventTime,
      #                                            g_list = g_list,
      #                                            t_list = t_list,
      #                                            N_g_list = N_g_list,
      #                                            use_last_treated_only = use_last_treated_only)
    } else if (estimand == "all"){

      A_0_list_simple <- create_A0_list_for_simple_average_ATE(g_list = g_list,
                                                               t_list = t_list,
                                                               N_g_list = N_g_list,
                                                               use_last_treated_only = use_last_treated_only)

      A_0_list_cohort <- create_A0_list_for_cohort_average_ATE(g_list = g_list,
                                                               t_list = t_list,
                                                               N_g_list = N_g_list,
                                                               use_last_treated_only = use_last_treated_only)

      A_0_list_calendar <- create_A0_list_for_calendar_average_ATE(g_list = g_list,
                                                                   t_list = t_list,
                                                                   N_g_list = N_g_list,
                                                                   use_last_treated_only = use_last_treated_only)
      # Event Study A_0_list
      if(length(eventTime)>1){
        n_eventTime = length(eventTime)
        aux_list = NULL
        for(j in 1:n_eventTime){
          aux_list_new <- create_A0_list_for_event_study(eventTime = eventTime[j],
                                                         g_list = g_list,
                                                         t_list = t_list,
                                                         N_g_list = N_g_list,
                                                         use_last_treated_only = use_last_treated_only)
          aux_list <- rbind(aux_list, aux_list_new)
        }

        A_0_list_event_study <- apply(aux_list, 2, rbind)

        A_0_list_event_study <-  purrr::pmap(.l = list(A_0_list_event_study),
                                             .f = function(A0){ return(do.call(rbind, A0)) }
        )


      } else {
        A_0_list_event_study <- create_A0_list_for_event_study(eventTime = eventTime,
                                                               g_list = g_list,
                                                               t_list = t_list,
                                                               N_g_list = N_g_list,
                                                               use_last_treated_only = use_last_treated_only)
      }

      my_list <- rbind(A_0_list_simple, A_0_list_cohort,  A_0_list_calendar, A_0_list_event_study)
      A0var_list <- apply(my_list, 2, rbind)

      A_0_list <-  purrr::pmap(.l = list(A0var_list),
                               .f = function(A0){ return(do.call(rbind, A0)) }
      )

      estim <- c("all_simple","all_cohort","all_calendar", paste0("all_ES",eventTime))


    }
  }

  balance_checks_Xhat <- compute_balance_test(Ybar_g_list = Ybar_g_list,
                                              A_0_list = A_0_list,
                                              S_g_list = S_g_list,
                                              N_g_list = N_g_list)


  se_Xhat <- sqrt(diag(as.matrix(balance_checks_Xhat$Xvar)))
  resultsDF <- data.frame(Xhat = balance_checks_Xhat$Xhat,
                          se_Xhat = se_Xhat,
                          t_test = abs(balance_checks_Xhat$Xhat/se_Xhat),
                          pvalue_t = 2*stats::pnorm(-abs(balance_checks_Xhat$Xhat/se_Xhat)),
                          Wald_test_Xhat = balance_checks_Xhat$Wald_test_Xhat,
                          pvalue_Wald = stats::pchisq(balance_checks_Xhat$Wald_test_Xhat,
                                                      df = length( balance_checks_Xhat$Xhat),
                                                      lower.tail = FALSE),
                          N = balance_checks_Xhat$N,
                          fisher_pval = NA,
                          fisher_supt_pval = NA,
                          num_fisher_permutations = NA,

                          estimand = estim
  )

  Xvar = NULL


  permutation_t_test <- function(
    df,
    A_0_list){


    g_level_summaries <- compute_g_level_summaries(df, is_balanced = TRUE)
    Ybar_g_list <- g_level_summaries$Ybar_g_List
    S_g_list <- g_level_summaries$S_g_List
    N_g_list <- g_level_summaries$N_g_List
    g_list <- g_level_summaries$g_list
    t_list <- g_level_summaries$t_list

    balance_checks_Xhat <- compute_balance_test(Ybar_g_list = Ybar_g_list,
                                                A_0_list = A_0_list,
                                                S_g_list = S_g_list,
                                                N_g_list = N_g_list)

    se_Xhat <- sqrt(diag(as.matrix(balance_checks_Xhat$Xvar)))
    t_test <- abs(balance_checks_Xhat$Xhat/se_Xhat)
    as.matrix(t_test, nrow = 1)
  }


  ## Do FRT, if specified
  permuteTreatment2 <- function(df, i_g_table, seed = NULL){
    #This function takes a data.frame with columns i and g, and permutes the values of g assigned to i
    # The input i_g_table has the unique combinations of (i,g) in df, and is calculated outside for speed improvements

    #Draw a random permutation of the elements of first_period_df
    if(!is.null(seed)) set.seed(seed)
    n = base::NROW(i_g_table)
    randIndex <-
      sample.int(n = n,
                 size = n,
                 replace = FALSE)

    #Replace first_period_df$g with a permuted version based on randIndex
    i_g_table$g <- i_g_table$g[randIndex]

    #Merge the new treatment assignments back with the original
    df$g <- NULL
    df <- dplyr::left_join(df,
                           i_g_table,
                           by = c("i"))

    return(as.data.frame(df))
  }


  if(compute_fisher == TRUE){

    #Find unique pairs of (i,g). This will be used for computing the permutations
    # i_g_table <- df %>%
    #              dplyr::filter(t == min(t)) %>%
    #              dplyr::select(i,g)

    i_g_table <- df_processed %>%
      dplyr::filter(t == min(t))
    i_g_table <- i_g_table[,c("i","g")]

    #Now, we compute the FRT for each seed, permuting treatment for each one
    #We catch any errors in the FRT simulations, and throw a warning if at least one has an error
    #(using the remaining draws to calculate frt)


    if(!is.null(seed)) set.seed(see)
    seed_frt <-  1:num_fisher_permutations * base::floor(stats::rexp(1, rate = 1/500))

    FRTResults_bal <-
      purrr::map(.x = seed_frt, #1:num_fisher_permutations,
                 .f = purrr::possibly(
                   .f = ~permutation_t_test(df = permuteTreatment2(df = df_processed,
                                                                   i_g_table = i_g_table,
                                                                   seed = .x),
                                            A_0_list = A_0_list) ,
                   otherwise = NULL)
      ) %>%
      purrr::discard(base::is.null)

    FRTResults_bal <- base::t(sapply(FRTResults_bal, base::rbind))
    FRTResults_bal <- as.data.frame(FRTResults_bal)

    successful_frt_draws <- base::NROW(FRTResults_bal)




    FRTResults2 <- FRTResults_bal

    if(successful_frt_draws < num_fisher_permutations){
      warning("There was an error in at least one of the FRT simulations. Removing the problematic draws.")
    }

    # Comput sup-t tests
    FRT_t <- as.numeric(base::apply(FRTResults_bal, 1, max))
    max_t <- max(resultsDF$t_test)
    resultsDF$fisher_supt_pval <- mean( max_t < FRT_t )

    # Need to compute p-value per row
    resultsDF$fisher_pval <- base::rowMeans(apply(FRTResults_bal, 1, '>', resultsDF$t_test))


    resultsDF$num_fisher_permutations <- successful_frt_draws
  }

  resultsDF <- as.data.frame(resultsDF)

  if(return_full_vcv){
    Xvar = balance_checks_Xhat$Xvar
  }

  list_results <- list(resultsDF = resultsDF,
                       Xvar = Xvar)

  return(list_results)





}
