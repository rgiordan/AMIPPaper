# You can optinally use DefineMacro here to set latex macros to use R quantities
# in the text.


# Define LaTeX macros that will let us automatically refer
# to simulation and model parameters.
DefineMacro <- function(macro_name, value, digits=3) {
  #sprintf_code <- paste("%0.", digits, "f", sep="")
  value_string <- format(value, big.mark=",", digits=digits, scientific=FALSE)
  cat("\\newcommand{\\", macro_name, "}{", value_string, "}\n", sep="")
}


DefineMacro("SimNumObs", sim_env$grid_list$num_obs, digits=0)
DefineMacro("SimTrueTheta", sim_env$grid_list$theta0, digits=1)

DefineMacro("SimAccNumObs", sim_env$acc_list$num_obs, digits=0)
DefineMacro("SimAccSigx", sim_env$acc_list$sig_x, digits=0)
DefineMacro("SimAccSigeps", sim_env$acc_list$sig_eps, digits=0)
DefineMacro("SimAccPercentMax", 100 * sim_env$acc_list$alpha_max, digits=0)
