



GetTreatmentInfluence <- function(regression_df, outcome, fmla,
                                  weights, lm_result) {
  # First, get the model matrices.
  y_all <- regression_df[[outcome]]
  keep_rows <- !is.na(y_all)
  y <- y_all[keep_rows]
  x <- as.array(model.matrix(fmla, regression_df))
  w <- weights[keep_rows]

  # Sanity checks that we're constructing the model matrices correctly.
  py_main$y <- y
  py_main$x <- as.array(x)
  py_main$w <- as.array(w)
  reticulate::py_run_string("
beta_py = regsens_rgiordandev.reg(y=y, x=x, w=w)
")
  stopifnot(max(abs(py_main$beta_py - lm_result$estimate)) < 1e-8)

  # A super redundant sanity check, but it's fast.
  beta <- solve(t(x * w) %*% x, t(x * w) %*% y)
  stopifnot(max(abs(beta - lm_result$estimate)) < 1e-8)

  # Get leverage.
  influence_matrix <- GetInfluence(y=y, x=x, w=w, beta=beta)
  treat_col <- which(names(regression_df) == "treatment")
  print(treat_col)
  influence_vec <- influence_matrix[treat_col, ]
  reg_list <- RegressionTargetChange(
      reg=lm_result, influence_vec=influence_vec, target_index=treat_col)

  return(list(influence_vec=influence_vec,
              reg_list=reg_list,
              keep_rows=keep_rows,
              treat_col=treat_col))
}



##### ITT
replicate_QJE_table_ITT <- function(outcome, dataset){
  # So that leverage rows index into a variable in the global scope, I
  # now do this filter outside the function.
  # dataset <- dataset %>%
  #   filter(sample_12m_resp == 1) # This conditional seemed to be everywhere in the do file.

  # Creating model formula
  fmla <- as.formula(paste0(outcome, "~ . + 0  - household_id"))
  # dplyr NSE
  # Only controls, treatment and outcome
  regression_df <- dataset %>%
    select_at(vars(treatment,
                   household_id, # Need these for clustering SEs
                   starts_with("ddd"), # These are the fixed effects
                   outcome))
  # Hardcoded weight according to Stata do file for survey data.
  # RJG: weight_12m is "Final 12-month weights" from raw_dta_surveys[[3]]
  model <- lm(data = regression_df,
              formula = fmla,
              weights = dataset$weight_12m)
  # Cluster on household
  clustered_SE <- vcovCL(model,
                         cluster = ~household_id) %>%
    diag() %>%
    sqrt()
  # Adding clustered standard errors, unadjusted p values and tidying.
  # Unlike the original, use the same names for statistics and std.error as tidy().
  tidy_model <- tidy(model,
                     quick = TRUE) %>% # quick means we compute our own SEs
    mutate(std.error = clustered_SE,
           statistic = estimate / std.error,
           pval = 2*pt(-abs(statistic), df = model$df))

  # Get influence
  infl_list <- GetTreatmentInfluence(
    regression_df, outcome, fmla, weights=dataset$weight_12m,
    lm_result=tidy_model)

  ExpandRows <- function(vec) {
    new_vec <- rep(NA, nrow(dataset))
    new_vec[infl_list$keep_rows] <- vec
    return(new_vec)
  }

  regression_output_df <- data.frame(
    influence=ExpandRows(infl_list$influence_vec),
    residuals=ExpandRows(model$residuals),
    yhat=ExpandRows(model$fitted.values),
    outcome=regression_df[[outcome]],
    treatment=regression_df$treatment)

  return(list(
    tidy_model=tidy_model,
    regression_df=regression_df,
    regression_output_df=regression_output_df,
    infl_list=infl_list,
    fmla=fmla
  ))
}
