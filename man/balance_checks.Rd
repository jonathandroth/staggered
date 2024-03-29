% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance_checks.R
\name{balance_checks}
\alias{balance_checks}
\title{Calculate balance checks in staggered rollout designs}
\usage{
balance_checks(
  df,
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
  seed = NULL
)
}
\arguments{
\item{df}{A data frame containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)}

\item{i}{The name of column containing the individual (cross-sectional unit) identifier. Default is "i".}

\item{t}{The name of the column containing the time periods. Default is "t".}

\item{g}{The name of the column containing the first period when a particular observation is treated, with Inf denoting never treated. Default is "g".}

\item{y}{The name of the column containing the outcome variable. Default is "y".}

\item{estimand}{The estimand to be calculated: "simple" averages all treated (t,g) combinations with weights proportional to N_g; "cohort" averages the ATEs for each cohort g, and then takes an N_g-weighted average across g; "calendar" averages ATEs for each time period, weighted by N_g for treated units, and then averages across time. "eventstudy" returns the average effect at the ''event-time'' given in the parameter EventTime.  The parameter can be left blank if a custom parameter is provided in A_theta_list. The argument is not case-sensitive.}

\item{A_0_list}{This parameter allow for specifying the matrices used to construct the Xhat vector of pre-treatment differences. If left NULL, the default is to use the scalar set of controls used in Callaway and Sant'Anna. If use_DiD_A0 = FALSE, then it uses the full vector possible comparisons of (g,g') in periods t<g,g'.}

\item{eventTime}{If using estimand = "eventstudy", specify what eventTime you want the event-study parameter for. The default is 0, the period in which treatment occurs. If a vector is provided, estimates are returned for all the event-times in the vector.}

\item{use_DiD_A0}{If this parameter is true, then Xhat corresponds with the scalar used by Callaway and Sant'Anna, so the Callaway and Sant'Anna estimator corresponds with beta=1. If it is false, the Xhat is a vector with all possible comparisons of pairs of cohorts before either is treated. The latter option should only be used when the number of possible comparisons is small relative to sample size.}

\item{use_last_treated_only}{If true, then A_0_list and A_theta_list are created to only make comparisons with the last treated cohorts (as suggested by Sun and Abraham), rather than using not-yet-treated units as comparisons. If set to TRUE (and use_DiD_A0 = TRUE), then beta=1 corresponds with the Sun and Abraham estimator.}

\item{compute_fisher}{If true, computes a Fisher Randomization Test using the studentized estimator.}

\item{num_fisher_permutations}{The number of permutations to use in the Fisher Randomization Test (if compute_fisher = TRUE). Default is 500.}

\item{return_full_vcv}{If this is true, then the function returns a list containing the full variance-covariance matrix for all Xhats.}

\item{skip_data_check}{If true, skips checks that the data is balanced and contains the colums i,t,g,y. Used in internal recursive calls to increase speed, but not recommended for end-user.}

\item{seed}{Set seed for permutations}
}
\value{
resultsDF A data.frame containing: estimate (the point estimate), se (the standard error), and se_neyman (the Neyman standard error). If a vector-valued eventTime is provided, the data.frame contains multiple rows for each eventTime and an eventTime column. If return_full_vcv = TRUE and estimand = "eventstudy", the function returns a list containing resultsDF and the full variance covariance for the event-study estimates (vcv) as well as the Neyman version of the covariance matrix (vcv_neyman). (If return_matrix_list = TRUE, it likewise returns a list containing lists of matrices used in the vcv calculation.)
}
\description{
This functions calculates Wald-tests for balance in staggered rollout designs proposed by Roth and Sant'Anna.
}
\examples{
set.seed(1234)
# load the officer data and subset it
df <- pj_officer_level_balanced
group_random <- sample(unique(df$assigned), 3)
df <- df[df$assigned \%in\% group_random,]
# Calculate balance checks for simple aggregation
balance_checks(df = df,
  i = "uid",
  t = "period",
  g = "first_trained",
  y = "complaints",
  estimand = "simple")
# Calculate balance checks for the cohort weighted average
balance_checks(df = df,
  i = "uid",
  t = "period",
  g = "first_trained",
  y = "complaints",
  estimand = "cohort")
# Calculate balance checks for the calendar weighted average
balance_checks(df = df,
  i = "uid",
  t = "period",
  g = "first_trained",
  y = "complaints",
  estimand = "calendar")
# Calculate balance checks for event-study aggregation for the first 24 months
# (month 0 is instantaneous effect)
event_bal_checks <- balance_checks(df = df,
  i = "uid",
  t = "period",
  g = "first_trained",
  y = "complaints",
  estimand = "eventstudy",
  eventTime = 0:23)
head(event_bal_checks)

}
\references{
\cite{Roth, Jonatahan, and Sant'Anna, Pedro H. C. (2021),
'Efficient Estimation for Staggered Rollout Designs', arXiv: 2102.01291, \url{https://arxiv.org/abs/2102.01291}.}
}
