# Main internal helper
create_A0_list_for_ATE_tg <- function(t,
                                      g,
                                      g_list,
                                      t_list,
                                      N_g_list,
                                      use_last_treated_only = FALSE){

  #if(t < g){warning("t is less than g. ATE(t,g) is zero by assumption")}
  if(t >= base::max(g_list)){base::stop("t is greater than max(g)-1; ATE(t,g) is not identified.")}

  # Create A0s for ATT_{t,g}
  # ------------------------

  treated_cohort_index <- base::which(g_list == g)
  if(!use_last_treated_only){
    #Create a list of which cohorts are eligible to be controls for each of the cohorts
    #This will be a null list if not eligible
    control_cohort_indices <- base::which(g_list > pmax(t,g))
  }else{
    #If use_last_treated_only, compare only to the last treated cohort (i.e. max(G))
    control_cohort_indices <- base::which((g_list > t) & (g_list ==  base::max(g_list)))
  }
  N_control_cohorts <- base::sum(N_g_list[control_cohort_indices])

  A_theta_t <- base::numeric(base::length(g_list))
  # weights for cohort in control_g_index
  A_theta_t[control_cohort_indices] <- - N_g_list[control_cohort_indices] / N_control_cohorts
  # Treated units
  A_theta_t[treated_cohort_index] <- A_theta_t[treated_cohort_index] + 1
  return(A_theta_t)
}


create_A0_list_for_event_study <- function(eventTime,
                                           g_list,
                                           t_list,
                                           N_g_DT,
                                           use_last_treated_only = FALSE){

  # Create A0s for an ``event-study'' coefficient at lag eventTime
  # This estimand is the average treatment effects for units eventTime periods after first being treated
  # The average is taken over all cohorts g such that there is some untreated cohort at g+eventTime
  # Cohorts are weighted by the cohort size (N_g)
  # This function returns the pre-period differences used as control in CS
  # i.e, the estimate of the the difference in period g-1 btwn cohort g and units not-yet-treated at g+eventTime
  # If use_last_treated_only = TRUE, then only the last cohort is used as a control

  maxG <- base::max(g_list)
  eligible_cohort_index <- base::which( ( pmax(g_list + eventTime, g_list) < maxG ) & ( pmax(g_list + eventTime, g_list) <= base::max(t_list) ) )
  N_g_list <- N_g_DT$N_g

  if(base::length(eligible_cohort_index) == 0){
    base::stop("There are no comparison cohorts for the given eventTime")
  }

  A_0 <- base::matrix(0, base::length(g_list), base::length(t_list))
  N_eligible <- sum(N_g_list[eligible_cohort_index])

  # length-G vector with weights for t, g, scaled by N_g/sum N_g
  create_A0_tg_helper <- function(gIndex) {
    g <- g_list[gIndex]
    (N_g_list[gIndex]/N_eligible) *
      create_A0_list_for_ATE_tg(g+eventTime, g, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # the results are placed in the time period columns corresponding to
  # the eligible cohorts offset by eventTime, minus 1
  tindices <- base::unlist(base::sapply(g_list[eligible_cohort_index], function(g) base::which(t_list == g-1)))

  # we can mechanically run the helper function for all eligible cohorts,
  # but not  all eligible cohorts will have a corresponding time column
  # (this is particularly plausible when the event time offset is large)
  gindices <- base::unlist(base::sapply(eligible_cohort_index, function(g) if (base::any(t_list == g_list[g]-1)) g else NULL))
  A_0[,tindices] <- base::sapply(eligible_cohort_index, create_A0_tg_helper)[,gindices]

  return(A_0)
}


create_A0_list_for_calendar_average_ATE <- function(g_list,
                                                    t_list,
                                                    N_g_DT,
                                                    use_last_treated_only = FALSE){

  t_eligible_index <- base::which(t_list >= base::min(g_list) & t_list < base::max(g_list))
  T_eligible <- base::length(t_eligible_index)
  N_g_list <- N_g_DT$N_g
  A_0 <- base::matrix(0, base::length(g_list), base::length(t_list))

  # G x T matrix; scaled internally
  create_A0_tg_helper <- function(t) {
      create_A0_list_for_ATE_calendar_t(t, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # t is fixed but the output is a matrix because across cohorts, the
  # relevant G x 1 vector is placed in the time column corresponding
  # to the g-1 for the gth cohort; hence create_A0_tg_helper returns a
  # matrix with those columns populated, for each t. We take the average
  # across time (add them and divide by T_eligible
  A_0 <- base::Reduce(`+`, base::lapply(t_list[t_eligible_index], create_A0_tg_helper))/T_eligible

  return(A_0)
}


create_A0_list_for_ATE_calendar_t <- function(t,
                                              g_list,
                                              t_list,
                                              N_g_list,
                                              use_last_treated_only = FALSE){

  treated_by_t_indices <- base::which(g_list <= t)
  N_total_treated <- base::sum(N_g_list[treated_by_t_indices])
  A_0_t <- base::matrix(0, base::length(g_list), base::length(t_list))

  # G x 1 vector containing weights for t, g; scaled by N_g / sum N_g
  create_A0_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_treated) *
          create_A0_list_for_ATE_tg(t, g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }

  # Note each vector is placed in the time column corresponding to
  # g-1, for some cohort g. The collection of vectors correspond to
  # all the eligible cohorts, but they have to be place in the time
  # columns corresponding to g-1.
  tindices <- base::unlist(base::sapply(g_list[treated_by_t_indices], function(g) base::which(t_list == g-1)))
  gindices <- base::unlist(base::sapply(treated_by_t_indices, function(g) if (base::any(t_list == g_list[g]-1)) g else NULL))
  A_0_t[,tindices] <- base::sapply(treated_by_t_indices, create_A0_tg_helper)[,gindices]

  return(A_0_t)
}


create_A0_list_for_cohort_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_DT,
                                                  use_last_treated_only = FALSE){

  N_g_list <- N_g_DT$N_g
  g_eligible_index <-  base::which(g_list < base::max(g_list) & g_list <= base::max(t_list))
  N_total_eligible <- base::sum(base::unlist(N_g_list[g_eligible_index]))

  # G x 1 vector; scaled by N_g / sum_g N_g
  create_A0_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_eligible) *
          create_A0_list_for_ATE_cohort_g(g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }

  # For each g, create_A0_list_for_ATE_cohort_g aggregates G x 1 vectors
  # across every t (takes the average). The resulting matrix is a
  # collection of vectors corresponding to each cohort, which must then
  # be placed into the time columns corresponding to each cohort g-1.
  A_0 <- base::matrix(0, base::length(g_list), base::length(t_list))
  tindices <- base::unlist(base::sapply(g_list[g_eligible_index], function(g) base::which(t_list == g-1)))
  gindices <- base::unlist(base::sapply(g_eligible_index, function(g) if (base::any(t_list == g_list[g]-1)) g else NULL))
  A_0[,tindices] <- base::sapply(g_eligible_index, create_A0_tg_helper)[,gindices]
  return(A_0)
}


create_A0_list_for_ATE_cohort_g <- function(g,
                                            g_list,
                                            t_list,
                                            N_g_list,
                                            use_last_treated_only = FALSE){
  treated_period_indices <- base::which((t_list >= g) & (t_list < base::max(g_list)))
  T_treated <- base::length(treated_period_indices)
  A_theta_g <- base::matrix(0, base::length(g_list), base::length(t_list))

  # G x 1 vector containing weights for t, g; will be scaled by parent function
  create_A0_tg_helper <- function(t) {
      create_A0_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }

  # For a given given G all vectors correspond to the same time column;
  # they will be placed in the g-1 time column by the parent function,
  # and here we only need to take the average.
  A_0_g <- base::Reduce(`+`, base::lapply(t_list[treated_period_indices], create_A0_tg_helper))/T_treated

  return(A_0_g)
}


create_A0_list_for_simple_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_DT,
                                                  use_last_treated_only = FALSE){

  N_g_list <- N_g_DT$N_g

  #Create a df with all the (g,t) pairs for which ATE is identified
  #Join in N_g for each of these pairs
  gt_df <- N_g_DT[data.table::CJ(g=g_list, t=t_list)[t >= g & t < base::max(g_list)]]

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- gt_df[,base::sum(N_g)]
  A_0 <- matrix(0, base::length(g_list), base::length(t_list))
  gt_df[,N_g := N_g / N_total]

  #sort by t
  data.table::setkey(gt_df, t, g)

  # generate the requisite vector for each t, g pair;
  # NB: This only returns a vector, and for a given t, g pair this is
  # the (t_list == (g-1))th column
  create_A0_tg_helper <- function(t, g, N_g) {
      N_g * create_A0_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }
  A_0_gt <- gt_df[, base::as.list(create_A0_tg_helper(t, g, N_g)), by=c("t", "g")]
  data.table::setnames(A_0_gt, c("t", "g", g_list))

  # sum all the vectors corresponding to the (t_list == (g-1))th column
  A_0_t    <- A_0_gt[,base::lapply(.SD, base::sum), by="g", .SDcols=base::as.character(g_list)]
  tindices <- base::unlist(base::sapply(A_0_t$g, function(g) base::which(t_list == g-1)))
  gindices <- base::unlist(base::sapply(A_0_t$g, function(g) if (base::any(t_list == g-1)) g else NULL))
  data.table::setkey(A_0_t, g)

  # place into G x T matrix
  A_0[,tindices] <- base::t(A_0_t[.(gindices), !"g"])
  return(A_0)
}
