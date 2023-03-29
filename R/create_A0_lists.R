# Main internal helper
create_A0_list_for_ATE_tg <- function(t,
                                      g,
                                      g_list,
                                      t_list,
                                      N_g_list,
                                      use_last_treated_only = FALSE){

  numPeriods <- length(t_list)
  #if(t < g){warning("t is less than g. ATE(t,g) is zero by assumption")}
  if(t >= max(g_list)){stop("t is greater than max(g)-1; ATE(t,g) is not identified.")}

  #Create A0s for ATT_{t,g}
  treated_cohort_index <- which(g_list == g)
  # event_time_index <- which(t_list == (g-1))

  #Create A0 for the treated units
  #NB: The old function looped through g_list and placed a 1 in the row
  #where (g == g_list) and column where (t == t_list). But you can just
  #put a one there directly.

  #NB: The old code made as many indices as number of cohorts, but the indices
  #were all identical, so there is no reason to make copies,
  if(!use_last_treated_only){
    #Create a list of which cohorts are eligible to be controls for each of the cohorts
    #This will be a null list if not eligible
    control_cohort_indices <- which(g_list > t)
  }else{
    #If use_last_treated_only, compare only to the last treated cohort (i.e. max(G))
    control_cohort_indices <- which((g_list > t) & (g_list ==  max(g_list)))
  }
  N_control_cohorts <- sum(N_g_list[control_cohort_indices])

  #This function creates the weights for the cohort in control_g_index
  #as a controls for the cohort in treated_g_index
  #
  #NB: The old function looped through all g indices and returned a row of 0s
  #if g is not in control_cohort_indices, else it returns a row of 0s with
  #N_g/N_control in the ((g-1) == t_list)th column. But, again, you can just put
  #that into those rows and column directly

  A_theta_t <- numeric(length(g_list))
  A_theta_t[control_cohort_indices] <- - N_g_list[control_cohort_indices] / N_control_cohorts
  A_theta_t[treated_cohort_index]   <- A_theta_t[treated_cohort_index] + 1
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

  maxG <- max(g_list)
  eligible_cohort_index <- which( ((g_list + eventTime) < maxG ) & ((g_list + eventTime) <= max(t_list) ) )
  N_g_list <- N_g_DT$N_g

  if(length(eligible_cohort_index) == 0){
    stop("There are no comparison cohorts for the given eventTime")
  }

  A_0 <- matrix(0, length(g_list), length(t_list))
  N_eligible <- sum(N_g_list[eligible_cohort_index])

  create_A0_tg_helper <- function(gIndex) {
    g <- g_list[gIndex]
    (N_g_list[gIndex]/N_eligible) *
      create_A0_list_for_ATE_tg(g+eventTime, g, g_list, t_list, N_g_list, use_last_treated_only)
  }
  tindices <- unlist(sapply(g_list[eligible_cohort_index], function(g) which(t_list == g-1)))
  gindices <- unlist(sapply(eligible_cohort_index, function(g) if (any(t_list == g_list[g]-1)) g else NULL))
  A_0[,tindices] <- sapply(eligible_cohort_index, create_A0_tg_helper)[,gindices]

  return(A_0)
}


create_A0_list_for_calendar_average_ATE <- function(g_list,
                                                    t_list,
                                                    N_g_DT,
                                                    use_last_treated_only = FALSE){

  t_eligible_index <- which(t_list >= min(g_list) & t_list < max(g_list))
  T_eligible <- length(t_eligible_index)
  N_g_list <- N_g_DT$N_g
  A_0 <- matrix(0, length(g_list), length(t_list))

  create_A0_tg_helper <- function(t) {
      create_A0_list_for_ATE_calendar_t(t, g_list, t_list, N_g_list, use_last_treated_only)
  }
  A_0 <- Reduce(`+`, lapply(t_list[t_eligible_index], create_A0_tg_helper))/T_eligible

  return(A_0)
}


create_A0_list_for_ATE_calendar_t <- function(t,
                                              g_list,
                                              t_list,
                                              N_g_list,
                                              use_last_treated_only = FALSE){

  treated_by_t_indices <- which(g_list <= t)
  N_total_treated <- sum(N_g_list[treated_by_t_indices])
  A_0_t <- matrix(0, length(g_list), length(t_list))

  create_A0_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_treated) *
          create_A0_list_for_ATE_tg(t, g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }
  tindices <- unlist(sapply(g_list[treated_by_t_indices], function(g) which(t_list == g-1)))
  gindices <- unlist(sapply(treated_by_t_indices, function(g) if (any(t_list == g_list[g]-1)) g else NULL))
  A_0_t[,tindices] <- sapply(treated_by_t_indices, create_A0_tg_helper)[,gindices]

  return(A_0_t)
}


create_A0_list_for_cohort_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_DT,
                                                  use_last_treated_only = FALSE){

  N_g_list <- N_g_DT$N_g
  g_eligible_index <-  which(g_list < max(g_list) & g_list <= max(t_list))
  N_total_eligible <- sum(unlist(N_g_list[g_eligible_index]))
  create_A0_tg_helper <- function(gIndex) {
      (N_g_list[gIndex] / N_total_eligible) *
          create_A0_list_for_ATE_cohort_g(g_list[gIndex], g_list, t_list, N_g_list, use_last_treated_only)
  }
  A_0 <- matrix(0, length(g_list), length(t_list))
  tindices <- unlist(sapply(g_list[g_eligible_index], function(g) which(t_list == g-1)))
  gindices <- unlist(sapply(g_eligible_index, function(g) if (any(t_list == g_list[g]-1)) g else NULL))
  A_0[,tindices] <- sapply(g_eligible_index, create_A0_tg_helper)[,gindices]
  return(A_0)
}


create_A0_list_for_ATE_cohort_g <- function(g,
                                            g_list,
                                            t_list,
                                            N_g_list,
                                            use_last_treated_only = FALSE){
  treated_period_indices <- which((t_list >= g) & (t_list < max(g_list)))
  T_treated <- length(treated_period_indices)
  A_theta_g <- matrix(0, length(g_list), length(t_list))
  create_A0_tg_helper <- function(t) {
      create_A0_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }
  A_0_g <- Reduce(`+`, lapply(t_list[treated_period_indices], create_A0_tg_helper))/T_treated
  return(A_0_g)
}


create_A0_list_for_simple_average_ATE <- function(g_list,
                                                  t_list,
                                                  N_g_DT,
                                                  use_last_treated_only = FALSE){

  N_g_list <- N_g_DT$N_g

  #Create a df with all the (g,t) pairs for which ATE is identified
  #Join in N_g for each of these pairs
  gt_df <- N_g_DT[data.table::CJ(g=g_list, t=t_list)[t >= g & t < max(g_list)]]

  #Calculate sum of N_g for all eligible (t,g) pairs
  N_total <- gt_df[,sum(N_g)]
  A_0 <- matrix(0, length(g_list), length(t_list))
  gt_df[,N_g := N_g / N_total]

  #sort by t
  data.table::setkey(gt_df, t, g)

  # generate the requisite vector for each t, g pair;
  # NB: This only returns a vector, and for a given t, g pair this is
  # the (t_list == (g-1))th column
  create_A0_tg_helper <- function(t, g, N_g) {
      N_g * create_A0_list_for_ATE_tg(t, g, g_list, t_list, N_g_list, use_last_treated_only)
  }
  A_0_gt <- gt_df[, as.list(create_A0_tg_helper(t, g, N_g)), by=c("t", "g")]
  data.table::setnames(A_0_gt, c("t", "g", g_list))

  # sum all the vectors corresponding to the (t_list == (g-1))th column
  A_0_t    <- A_0_gt[,lapply(.SD, sum), by="g", .SDcols=as.character(g_list)]
  tindices <- unlist(sapply(A_0_t$g, function(g) which(t_list == g-1)))
  gindices <- unlist(sapply(A_0_t$g, function(g) if (any(t_list == g-1)) g else NULL))
  data.table::setkey(A_0_t, g)

  # place into G x T matrix
  A_0[,tindices] <- t(A_0_t[.(gindices), !"g"])
  return(A_0)
}
