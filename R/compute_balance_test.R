#' @title Wald-test if Xhat is statistically different from zero
#' @description \code{balance_check} computes the Wald test-statistic (F-test) for the null that E\[X\]=0
#' @param Ybar_g_list Ybar_g_list
#' @param A_0_list A_0_list
#' @param S_g_list S_g_list
#' @param N_g_list N_g_list
#' @return Wald_test_Xhat Wald-test statistic for the balance test.
compute_balance_test <- function(Ybar_g_list,
                                 A_0_list,
                                 S_g_list,
                                 N_g_list
){

  # First compute Xhat
  Xhat <- compute_Xhat(Ybar_g_list,
                       A_0_list)
  # Compute variance of Xhat
  Xvar_list <- purrr::pmap(.l = list(A_0_list,
                                     S_g_list, N_g_list),
                           # .f = function(A0,S,N){ return(1/N * A0 %*% S %*% base::t(A0) )  }
                           .f = function(A0,S,N){ return(1/N * eigenMapMatMult( eigenMapMatMult(A0,S) , base::t(A0) ) ) }
  )
  Xvar <- base::Reduce(f = '+',
                       x= Xvar_list)

  Xvar <- base::as.matrix(Xvar)
  # Get the inverse
  Xvar_inv <- MASS::ginv(Xvar)
  # Compute sample size
  sample_size <- purrr::reduce(N_g_list, base::sum)

  # Balance check
  Wald_test_Xhat <- base::t(Xhat) %*% Xvar_inv %*% Xhat


  return(list(Xhat = Xhat,
              Xvar = Xvar,
              N = sample_size,
              Wald_test_Xhat = Wald_test_Xhat
  )
  )
}
