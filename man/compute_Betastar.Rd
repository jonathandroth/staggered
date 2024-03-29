% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compute_efficient_estimator_and_se.R
\name{compute_Betastar}
\alias{compute_Betastar}
\title{Plug-in efficient Beta hat}
\usage{
compute_Betastar(
  Ybar_g_list,
  A_theta_list,
  A_0_list,
  S_g_list,
  N_g_list,
  Xvar_list = NULL
)
}
\arguments{
\item{Ybar_g_list}{Ybar_g_list}

\item{A_theta_list}{A_theta_list}

\item{A_0_list}{A_0_list}

\item{S_g_list}{S_g_list}

\item{N_g_list}{N_g_list}

\item{Xvar_list}{Xvar_list}
}
\value{
betastar Vector of plug-in efficient betahat estimates.
}
\description{
\code{compute_Betastar} computes the plug-in efficient betahat
}
