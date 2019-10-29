get_med_per_row <- function(prob_below_row, ks){

  if (prob_below_row[1] > 0.5) return(0)
  if (rev(prob_below_row)[1] < 0.5) return(Inf)

  f_row <- approxfun(log(ks), qlogis(unlist(prob_below_row)))
  exp(uniroot(f_row, lower = log(ks[1]), upper = log(rev(ks)[1]))$root)

}

##########################

get_med_for_given_col <- function(xi_col, eta_col, ks, sigma_error, ui, vi, samples = FALSE){

  pyks <- purrr::map_dfc(ks,
                         get_marg_pyk,
                         xi_col = xi_col,
                         eta_col = eta_col,
                         sigma_error = sigma_error,
                         ui = ui,
                         vi = vi,
                         samples = TRUE)

  med <- apply(pyks, 1, get_med_per_row, ks = ks)

  if (samples)  return(med)

  as.data.frame(t(quantile(med, c(0.025, 0.5, 0.975))))

}

##########################

marg_med <- function(newdata, nsims, ks, fit){

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
                        get_med_for_given_col,
                        ks = ks,
                        sigma_error = sigma_error,
                        ui = ui_vi$ui,
                        vi = ui_vi$vi))
}
