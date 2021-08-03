#' @title Procedural Justice Training Program in the Chicago Police Department
#'
#' @description Data from a large-scale procedural justice training program in the Chicago Police Department analyzed by
#' Wood, Tyler, Papachristos, Roth and Sant'Anna (2020) and Roth and Sant'Anna (2021). The data contains a balanced panel
#' of 7,785 police officers in Chicago who were randomly given a procedural justice training on different dates, and who
#'  remained in the police force throughout the study period (from January 2011 to December 2016).
#'
#' @format A data frame with 560520 observations (7,785 police officers and 72 months) and 12 variables:
#' \describe{
#'   \item{uid}{identifier for the police officer}
#'   \item{month}{month and year of the observation}
#'   \item{assigned}{month-year of first training assignment}
#'   \item{appointed}{appointment date}
#'   \item{resigned}{Date the police officer resigned. NA if he/she did not resigned by the time data was collected}
#'   \item{birth_year}{Officer's year of birth}
#'   \item{assigned_exact}{Exact date of first training assignment}
#'   \item{complaints}{Number of complaints (setlled and sustained)}
#'   \item{sustained}{Number of sustained complaints}
#'   \item{force}{Number of times force was used}
#'   \item{period}{Time period: 1 - 72}
#'   \item{first_trained}{Time period first exposed to treatment (Treatment cohort/group)}
#' }
#' @source Wood, Tyler, Papachristos, Roth and Sant'Anna (2020) and Roth and Sant'Anna (2021).
#' @references
#'   \cite{Roth, Jonatahan, and Sant'Anna, Pedro H. C. (2021),
#'   'Efficient Estimation for Staggered Rollout Designs', arXiv: 2102.01291, \url{https://arxiv.org/abs/2102.01291}.}
#'
#'
#'   \cite{Wood, George, Tyler, Tom R., Papachristos, Andrew P., Roth, Jonathan and Sant'Anna, Pedro H. C. (2020),
#'    'Revised findings for "Procedural justice training reduces police use of force and complaints against officers",
#'     \doi{10.31235/osf.io/xf32m}. }
"pj_officer_level_balanced"
