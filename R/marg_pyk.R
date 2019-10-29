get_marg_pyk <- function(xi_col, eta_col, k, sigma_error, ui, vi, samples = FALSE){
  pi <- plogis(c(xi_col) + ui)
  if (k == 0) {
    marg_pyk <- rowMeans(pi)
  }
  else {
    mu <- c(eta_col) + vi
    marg_pyk <- rowMeans(pi + (1 - pi) * pnorm((log(k) - mu) / sigma_error))
  }
  if (samples) return(marg_pyk)
  return(as.data.frame(t(quantile(marg_pyk, probs = c(0.025, 0.5, 0.975)))))
}


marg_pyk_q <- function(k, newdata, nsims, fit){

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
                        get_marg_pyk,
                        k = k,
                        sigma_error = sigma_error,
                        ui = ui_vi$ui,
                        vi = ui_vi$vi))
}
