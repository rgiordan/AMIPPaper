
import paragami
import vittles
import autograd
from autograd import numpy as np
from autograd import scipy as sp

#################################
# ADVI

def get_advi_pattern(param_dim):
    advi_param_pattern = paragami.PatternDict(free_default=True)
    advi_param_pattern['mean'] = paragami.NumericVectorPattern(param_dim)
    advi_param_pattern['log_sd'] = paragami.NumericVectorPattern(param_dim)
    return advi_param_pattern


def encode_draws(advi_params, base_draws):
    """Return a matrix of encoded draws, each column of which is a
    free parameter."""
    return \
        base_draws * np.exp(advi_params['log_sd'][:, None]) + \
        advi_params['mean'][:, None]


def get_advi_entropy(advi_params):
    return np.sum(advi_params['log_sd'])


def get_advi_loss(advi_params, base_draws, param_pattern, eval_log_prob):
    """eval_log_prob must be a function returning the target log density in the
    _contrained_ space.  The log abs jacobian is added here.
    """
    param_free_mat = encode_draws(advi_params, base_draws)
    num_draws = base_draws.shape[1]
    lp_sum = 0.0
    for d in range(num_draws):
        params = param_pattern.fold(param_free_mat[:, d], free=True)
        log_det_jac = param_pattern.log_abs_det_unfreeing_jacobian(params)
        lp_sum += eval_log_prob(params) + log_det_jac

    advi_entropy = get_advi_entropy(advi_params)
    return -1 * (lp_sum / num_draws + advi_entropy)


def get_base_draws(num_draws, param_pattern):
    param_dim = param_pattern.flat_length(free=True)
    return np.random.normal(size=(param_dim, num_draws))


def get_advi_objective(advi_param_pattern, base_draws, param_pattern, eval_log_prob):
    get_log_loss_free = paragami.FlattenFunctionInput(
        lambda advi_params: \
            get_advi_loss(advi_params, base_draws, param_pattern, eval_log_prob),
        patterns=advi_param_pattern,
        free=True)
    objective = paragami.OptimizationObjective(get_log_loss_free)
    return objective


def get_flat_draws(advi_params, base_draws, param_pattern):
    """Return a  matrix of flat draws which can be compared with MCMC output.
    """
    flat_advi_draws = []
    params_free_mat = encode_draws(advi_params, base_draws)
    for d in range(params_free_mat.shape[1]):
        params_flat = param_pattern.flatten(
            param_pattern.fold(params_free_mat[:, d], free=True), free=False)
        flat_advi_draws.append(params_flat)
    return np.array(flat_advi_draws)


def get_advi_mean(advi_params, base_draws, param_pattern):
    flat_draws = get_flat_draws(advi_params, base_draws, param_pattern)
    return np.mean(flat_draws, axis=0)
