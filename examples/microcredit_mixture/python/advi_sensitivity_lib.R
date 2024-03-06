
# Convert a list of parameter estimates into a tidy dataframe.
GatherMicrocreditResult <- function(result, param_names=NULL, increment_indices=FALSE) {
    # real mu[2];
    # real tau[2];
    # real<lower=0> sd_mu[2];
    # real<lower=0> sd_tau[2];
    # real sigma_control[2];
    # real sigma_TE[2];
    # real<lower=0> sd_sigma_control[2];
    # real<lower=0> sd_sigma_TE[2];
    # matrix[K,2] mu_k;
    # matrix[K,2] tau_k;
    # matrix[K,2] sigma_control_k;
    # matrix[K,2] sigma_TE_k;
    # matrix[M-1,P] beta; // the parent parameters minus the Mth category
    # matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
    # matrix[M,P] beta_k_raw[K]; // the hierarchical increments

    if (!is.matrix(result)) {
        result <- matrix(result, nrow=1)
    }
    if (!is.null(param_names)) {
        colnames(result) <- param_names
    }
    gathered_df <- tidybayes::gather_draws(
            model=result,
            mu[s], tau[s],
            sd_mu[s], sd_tau[s],
            sigma_control[s], sigma_TE[s],
            sd_sigma_control[s], sd_sigma_TE[s],
            mu_k[k, s], tau_k[k, s],
            sigma_control_k[k, s], sigma_TE_k[k, s],
            beta[m, p], sigma[m, p],
            beta_k_raw[k,m,p]) %>%
        select(-.chain, -.iteration) %>%
        rename(value=.value, row=.draw, param=.variable) %>%
        ungroup()
    if (increment_indices) {
        # Change from zero to one-based indexing.
        gathered_df <- gathered_df %>%
            mutate(k=k + 1, m=m + 1, p=p + 1, s=s + 1)
    }
    return(gathered_df)
}



#' 
#' # Convert the fit data into a ModelGrads object for zaminfluence.
#' #' @param mean_params:  A vector of the estimates of the main model parameters
#' #' @param mean_influence_mat:  The influence matrix for the main model parameters
#' #' @param lr_cov: The LR covariance matrix of the main model parameters
#' #' @param custom_advi_mean: A vector of estimates of functions of the main parameters
#' #' @param custom_influence_mat:The influence matrix for functions of the main parameters
#' #' @param custom_advi_mean: The LR covariance matrix for functions of the main parameters
#' #' 
#' #' @return a ModelGrads object for use with zaminfluence
#' PythonToModelGrads <- function(param_names,
#'                                mean_params,
#'                                mean_influence_mat,
#'                                lr_cov,
#'                                custom_param_names,
#'                                custom_advi_mean,
#'                                custom_influence_mat,
#'                                custom_advi_cov) {
#'     
#'     # Sanity check dimensions
#'     stopifnot(length(mean_params) == ncol(mean_influence_mat))
#'     stopifnot(length(custom_advi_mean) == ncol(custom_influence_mat))
#'     stopifnot(nrow(custom_influence_mat) == nrow(custom_influence_mat))
#'     stopifnot(ncol(mean_influence_mat) == ncol(lr_cov))
#'     stopifnot(ncol(custom_influence_mat) == ncol(custom_advi_cov))
#'     
#'     n_obs <- nrow(mean_influence_mat)
#'     
#'     param_names <- c(param_names, custom_param_names)
#'     betahat <- c(mean_params, custom_advi_mean)
#'     ses <- c(sqrt(diag(lr_cov)), sqrt(diag(custom_advi_cov)))
#'     infl_mat <- cbind(mean_influence_mat, custom_influence_mat)
#'     
#'     return(list(model_fit=NA,
#'                 n_obs=n_obs,
#'                 parameter_names=param_names,
#'                 
#'                 betahat=betahat,
#'                 se=ses,
#'                 weights=rep(1, n_obs),
#'                 se_group=NULL,
#'                 
#'                 beta_grad=t(infl_mat),
#'                 se_grad=matrix(0, nrow=ncol(infl_mat), ncol=n_obs),
#'                 
#'                 RerunFun=function(...) { stop("Not implemented") }
#'     ))
#' }
#' 


# Convert the fit data into a ModelGrads object for zaminfluence.
#' @param mean_params:  A vector of the estimates of the main model parameters
#' @param mean_influence_mat:  The influence matrix for the main model parameters
#' @param lr_cov: The LR covariance matrix of the main model parameters
#' @param custom_advi_mean: A vector of estimates of functions of the main parameters
#' @param custom_influence_mat:The influence matrix for functions of the main parameters
#' @param custom_advi_mean: The LR covariance matrix for functions of the main parameters
#' 
#' @return a ModelGrads object for use with zaminfluence
PythonToModelGrads <- function(param_names,
                               mean_params,
                               mean_influence_mat,
                               lr_cov,
                               custom_param_names,
                               custom_advi_mean,
                               custom_influence_mat,
                               custom_advi_cov) {
  
  # Sanity check dimensions
  stopifnot(length(param_names) == length(mean_params))
  stopifnot(ncol(lr_cov) == length(mean_params))
  stopifnot(ncol(mean_influence_mat) == length(mean_params))
  stopifnot(length(custom_advi_mean) == ncol(custom_influence_mat))
  stopifnot(nrow(custom_influence_mat) == nrow(custom_influence_mat))
  stopifnot(ncol(mean_influence_mat) == ncol(lr_cov))
  stopifnot(ncol(custom_influence_mat) == ncol(custom_advi_cov))
  
  n_obs <- nrow(mean_influence_mat)
  all_param_names <- c(param_names, custom_param_names)
  betahat <- c(mean_params, custom_advi_mean)
  ses <- c(sqrt(diag(lr_cov)), sqrt((custom_advi_cov)))
  infl_mat <- cbind(mean_influence_mat, custom_influence_mat)
  colnames(infl_mat) <- all_param_names
  
  model_fit <-
    ModelFit(fit_object=NA,
             num_obs=n_obs,
             param=betahat,
             se=ses,
             parameter_names=all_param_names,
             weights=NULL, se_group=NULL)
  
  se_grad <- t(infl_mat)
  se_grad[,] <- 0

  model_grads <- ModelGrads(
    model_fit,
    param_grad=t(infl_mat),
    se_grad=se_grad,
    RerunFun=function(...) { stop("Not implemented") })
  
  return(model_grads)
}



InitializePythonAndLoadFit <- function(initial_fit_filename, sens_filename, base_dir, python_path="python") {
    # Initialize Python
    reticulate::use_python(python_path)
    reticulate::py_run_string("
import autograd
import autograd.numpy as np
from copy import deepcopy
import paragami
import microcredit_mixture_rgiordandev
from microcredit_mixture_rgiordandev import microcredit_vb_lib as mm_lib
")
    py_main <- reticulate::import_main()
    py_main$initial_fit_filename <- initial_fit_filename
    
    # Load the ADVI results
    reticulate::py_run_string("
advi_dict = np.load(initial_fit_filename)
advi_param_pattern = paragami.get_pattern_from_json(str(advi_dict['advi_param_pattern_json']))
param_pattern = paragami.get_pattern_from_json(str(advi_dict['param_pattern_json']))
advi_params_free = advi_dict['advi_params_free']
advi_params = advi_param_pattern.fold(advi_dict['advi_params_free'], free=True)
base_draws = advi_dict['base_draws']
mean_params = mm_lib.get_advi_mean(advi_params, base_draws, param_pattern)
flat_names = param_pattern.flat_names(free=False, delim='')
hess0 = advi_dict['hess0']
advi_time = advi_dict['opt_time']
")
    
    # Load the sensitivity analysis
    py_main$sens_filename <- sens_filename
    reticulate::py_run_string("
sens_dict = np.load(sens_filename)
influence_mat = sens_dict['influence_mat']
mean_influence_mat = sens_dict['mean_influence_mat']
lr_cov = sens_dict['lr_cov']
sens_use_smuber = sens_dict['use_smuber']
")
    return(py_main)
}


# Compute custom functions of the advi parameters
DefineCustomParameterFunctions <- function(py_main) {
    # We assume that InitializePythonAndLoadFit() has been run.
    py_run_string("
def get_custom_param_vector(params):
    return np.log(params['sd_tau'])
")
    
    py_run_string("
def get_custom_advi_mean(advi_params):
    params_free_mat = mm_lib.encode_draws(advi_params, base_draws)
    advi_draws = [
        param_pattern.fold(params_free_mat[:, d], free=True) for
            d in range(params_free_mat.shape[1])
    ]
    param_draws = np.array([ 
        get_custom_param_vector(advi_draw) for advi_draw in advi_draws ])
    return np.mean(param_draws, axis=0)
")

    py_run_string("
def get_custom_advi_sensitivity(advi_params):
    custom_advi_mean = get_custom_advi_mean(advi_params)

    get_custom_advi_mean_free = paragami.FlattenFunctionInput(
        get_custom_advi_mean, advi_param_pattern, free=True, argnums=[0])
    custom_advi_mean_free_jacobian = autograd.jacobian(
        get_custom_advi_mean_free, argnum=0)(advi_params_free)
    custom_influence_mat = custom_advi_mean_free_jacobian @ influence_mat

    # Linear response covariance for the custom parameters.
    custom_lr_cov = \\
        custom_advi_mean_free_jacobian @ \\
        np.linalg.solve(hess0, custom_advi_mean_free_jacobian.T)

    return custom_advi_mean, custom_lr_cov, custom_influence_mat
")
}

EvalCustomParameters <- function(py_main, advi_params) {
    # Assume that DefineCustomParameterFunctions has been run
    py_main$this_advi_params <- advi_params
    py_run_string("
custom_advi_mean, custom_lr_cov, custom_influence_mat = \\
    get_custom_advi_sensitivity(this_advi_params)
")
    return(list(
        param_names=sprintf("log_sd_tau[%d]", 1:length(py_main$custom_advi_mean) - 1),
        mean=py_main$custom_advi_mean,
        sd=sqrt(diag(py_main$custom_lr_cov)),
        infl=py_main$custom_influence_mat))
}



# 
# GetTargetADVIGrads <- function(advi_influence, param_name,
#                                sig_num_ses=qnorm(0.975)) {
#     target_index <- which(advi_influence$param_names == param_name)
#     stopifnot(all(advi_influence$param_names == rownames(advi_influence$beta_grad)))
#     if (length(target_index) != 1) {
#         stop("Error finding target regressor in the regression.")
#     }
# 
#     grad_df <-
#         data.frame(
#             row=1:advi_influence$n_obs,
#             obs_per_row=1,
#             weights=advi_influence$weights,
#             se_grad=advi_influence$weights * advi_influence$se_grad[target_index,],
#             beta_grad=advi_influence$weights * advi_influence$beta_grad[target_index, ])
# 
#     attr(grad_df, "n_obs") <- advi_influence$n_obs
#     attr(grad_df, "n_grad_rows") <- advi_influence$n_obs
#     attr(grad_df, "data_row_cols") <- "row"
#     attr(grad_df, "obs_per_row_col") <- "obs_per_row"
#     base_vals <- c(
#         advi_influence$betahat[target_index],
#         advi_influence$se[target_index])
#     names(base_vals) <- c("beta", "se")
#     attr(grad_df, "base_vals") <- base_vals
#     attr(grad_df, "target_regressor") <- param_name
#     attr(grad_df, "target_index") <- target_index
# 
#     grad_df <- zaminfluence::AppendIntervalColumns(grad_df, sig_num_ses=sig_num_ses)
# 
#     return(grad_df)
# }

# 
# ProcessADVIInfluence <- function(mean_influence_mat,
#                                  param_name,
#                                  param_names,
#                                  base_means,
#                                  base_ses) {
#     n_obs <- nrow(mean_influence_mat)
#     # This is the expected output format from ComputeRegressionInfluence.
#     advi_influence <-
#         list(lm_result=NA,
#              n_obs=n_obs,
# 
#              param_names=param_names,
#              betahat=base_means,
#              se=base_ses,
#              weights=rep(1, n_obs),
# 
#              beta_grad=t(mean_influence_mat),
#              se_grad=matrix(0, nrow=ncol(mean_influence_mat), ncol=n_obs))
#     influence_dfs <-
#         GetTargetADVIGrads(advi_influence, param_name) %>%
#         SortAndAccumulate(advi_grads)
#     return(influence_dfs)
# }



