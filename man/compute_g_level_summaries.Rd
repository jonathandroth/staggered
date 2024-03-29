% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_efficient_estimator_and_se.R
\name{compute_g_level_summaries}
\alias{compute_g_level_summaries}
\title{Calculate group level summary statistics}
\usage{
compute_g_level_summaries(df, is_balanced = TRUE)
}
\arguments{
\item{df}{A data table containing panel data with the variables y (an outcome), i (an individual identifier), t (the period in which the outcome is observe), g (the period in which i is first treated, with Inf denoting never treated)}

\item{is_balanced}{If true, the df has previously been balanced so this does not need to be done internally.}
}
\value{
Y_bar_list A list of the means of the outcomes for each cohort g

S_g_list A list of covariance matrices for the outcomes for each cohort g

N_g_DT A data table of the number of observations for each cohort g

g_list A list of when the cohorts were first treated

t_list A list of the the time periods for the outcome. The vector of outcomes corresponds with this order.
}
\description{
This function computes the mean-vector and covariance matrix of the outcomes for each cohort, where a cohort g is a group of units first treated in period g
}
