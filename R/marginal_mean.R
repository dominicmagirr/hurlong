#' Simulate correlated random effects for a new patient
#'
#' Simulate matrices `ui` and `vi` of correlated random effects
#' from longitudinal hurdle model, where `ui` corresponds to the
#' zero part of the model, and `vi` corresponds to the non-zero part.
#' The number of rows is equal to the number of MCMC samples, the
#' number of columns is equal to the number of simulations used
#' in the inner Monte Carlo integration
#'
#' @param nsims The number of simulations used in the inner Monte Carlo integration.
#' @param sigma_u MCMC draws from the parameter `sigma_u`
#' @param sigma_v MCMC draws from the parameter `sigma_v`
#' @param rho MCMC draws from the parameter `rho`
#' @return A list of two matrices, `ui` and `vi`.
#' @export
#'
get_ui_vi <- function(nsims, sigma_u, sigma_v, rho){

  n_stan_sims <- length(sigma_u)

  ui <- matrix(rnorm(n_stan_sims * nsims, 0, sigma_u),
               nrow = n_stan_sims,
               ncol = nsims)

  cond_var_vi <- (1 - rho ^ 2) * sigma_v ^ 2
  cond_mean_vi <- ui * sigma_v / sigma_u * rho
  vi <- matrix(rnorm(n_stan_sims * nsims,
                     mean = cond_mean_vi,
                     sd = sqrt(cond_var_vi)),
               nrow = n_stan_sims,
               ncol = nsims)

  list(ui = ui, vi = vi)

}

#' Apply a function to MCMC draws of model parameters to give a sample from 'marginal mean'
#'
#' @param xi_col A vector of MCMC draws from the parameter xi
#' @param eta_col A vector of MCMC draws from the parameter eta
#' @param sigma_error MCMC draws from the parameter `sigma_error`
#' @param ui A matrix of draws from random effect (zero part)
#' @param vi A matrix of draws from random effect (non-zero part)
#' @param samples Return the samples rather than summarizing (default is FALSE)
#' @return Quantiles 0.025, 0.5, 0.975 from a posterior sample of marginal mean
#' @export
#'
get_marg_mean <- function(xi_col, eta_col, sigma_error, ui, vi, samples = FALSE){
  pi <- plogis(c(xi_col) + ui)
  mu <- c(eta_col) + vi
  marg_mean <- rowMeans((1 - pi) * exp(mu + sigma_error ^ 2  / 2))
  if (samples) return(marg_mean)
  return(as.data.frame(t(quantile(marg_mean, probs = c(0.025, 0.5, 0.975)))))

}


#' Summarize the posterior distribution of the 'marginal mean' at a give `newdata`
#'
#' @param newdata A data frame with covariates.
#' @param nsims The number of Monte Carlo samples used in the inner integration when calculating the marginal mean
#' @param fit The brmsfit object
#' @return Quantiles 0.025, 0.5, 0.975 from a posterior sample of marginal mean, evaluated at `newdata`
#' @export
#'
marg_mean_q <- function(newdata, nsims, fit){

  ps <- posterior_samples(fit)
  sigma_error <- ps[,"sigma"]
  sigma_u <- ps[,"sd_id__hu_Intercept"]
  sigma_v <- ps[,"sd_id__Intercept"]
  rho <- ps[,"cor_id__Intercept__hu_Intercept"]

  ui_vi <- get_ui_vi(nsims = nsims,
                     sigma_u = sigma_u,
                     sigma_v = sigma_v,
                     rho = rho)

  eta <- fitted(fit,
                newdata = newdata,
                re_formula = NA,
                dpar = "mu",
                summary = FALSE)

  xi <- qlogis(fitted(fit,
                      newdata = newdata,
                      re_formula = NA,
                      dpar = "hu",
                      summary = FALSE))

  cbind(newdata,
        purrr::map2_dfr(as.list(as.data.frame(xi)),
                        as.list(as.data.frame(eta)),
                        get_marg_mean,
                        sigma_error = sigma_error,
                        ui = ui_vi$ui,
                        vi = ui_vi$vi))
}






