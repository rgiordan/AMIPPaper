
import paragami
import vittles
import autograd
from autograd import numpy as np
from autograd import scipy as sp

def get_param_pattern(num_sites):
    # Priors:
    # mu, tau, sd_mu, sd_tau, sigma_k ~ P(.)
    #
    # Hierarchy:
    # mu_k | mu, sd_mu ~ P((mu_k - mu) / sd_mu)
    # tau_k | tau, sd_tau ~ P((tau_k - tau) / sd_tau)
    #
    # Data:
    # y_nk | ... ~ P((y_nk - (mu_k + treat * tau_k)) / sigma_k)
    num_sites = int(num_sites)
    param_pattern = paragami.PatternDict(free_default=True)

    param_pattern['mu'] = paragami.NumericScalarPattern()
    param_pattern['tau'] = paragami.NumericScalarPattern()
    param_pattern['sd_mu'] = paragami.NumericScalarPattern(lb=0.0)
    param_pattern['sd_tau'] = paragami.NumericScalarPattern(lb=0.0)

    param_pattern['mu_k'] = paragami.NumericVectorPattern(num_sites)
    param_pattern['tau_k'] = paragami.NumericVectorPattern(num_sites)
    param_pattern['sigma_k'] = paragami.NumericVectorPattern(num_sites, lb=0.0)

    return param_pattern


def get_prior_param_pattern():
    prior_pattern = paragami.PatternDict()
    def append_loc_scale(param_name):
        prior_pattern[param_name + '_loc'] = paragami.NumericScalarPattern()
        prior_pattern[param_name + '_scale'] = \
            paragami.NumericScalarPattern(lb=0.0)

    append_loc_scale('mu')
    append_loc_scale('tau')
    append_loc_scale('sd_mu')
    append_loc_scale('sd_tau')
    append_loc_scale('sigma_k')

    return prior_pattern


def get_default_prior_params():
    prior_params = dict()
    def append_loc_scale(param_name, loc, scale):
        prior_params[param_name + '_loc'] = loc
        prior_params[param_name + '_scale'] = scale

    append_loc_scale('mu', 0, 100)
    append_loc_scale('tau', 0, 100)
    append_loc_scale('sd_mu', 0, 100)
    append_loc_scale('sd_tau', 0, 100)
    append_loc_scale('sigma_k', 0, 100)
    return prior_params

# densitites

def get_normal_lpdf(z, mu, sigma):
    info = (1 / sigma**2)
    lp = -0.5 * info * (z - mu)**2 + 0.5 * np.log(info)
    # print('lp shape:', lp.shape,
    #       'z shape:', z.shape,
    #       'mu shape:', mu.shape,
    #       'sigma shape:', sigma.shape)
    return lp

def get_cauchy_lpdf(z, loc, scale):
    z_norm = (z - loc) / scale
    return -np.log(scale) - np.log(1 + z_norm ** 2)

# Smuber = "Smooth Huber".  The -log(sigma) comes from the rescaling Jacobian.
def get_smuber_lpdf(z, mu, sigma):
    z_norm = (z - mu) / sigma
    z_norm2 = z_norm ** 2
    return -1 * z_norm2 / np.sqrt(1 + z_norm2) - np.log(sigma)


def get_regression_lpdf(xtx, xty, yty, n_obs, beta, y_info):
    """Regression log likelihood for y_n ~ N(x_n^T beta, 1 / y_info)
    in terms of sufficient statistics.
    """

    x_dim = xtx.shape[0]
    assert xtx.shape[0] == xtx.shape[1]
    assert x_dim == len(beta)
    assert x_dim == len(xty)

    return \
        -0.5 * y_info * (
            yty -
            2 * np.dot(beta, xty) +
            np.trace(xtx @ np.outer(beta, beta))) + \
        0.5 * n_obs * np.log(y_info)


# Modeling functions

def get_log_prior(params, prior_params, config):
    def get_normal_prior(param_name):
        return np.sum(get_normal_lpdf(
            params[param_name],
            prior_params[param_name + '_loc'],
            prior_params[param_name + '_scale']))

    def get_cauchy_prior(param_name):
        return np.sum(get_cauchy_lpdf(
            params[param_name],
            prior_params[param_name + '_loc'],
            prior_params[param_name + '_scale']))

    lp = 0.0
    lp += get_normal_prior('mu')
    lp += get_normal_prior('tau')
    lp += get_cauchy_prior('sd_mu')
    lp += get_cauchy_prior('sd_tau')
    lp += get_cauchy_prior('sigma_k')

    return lp


def get_hierarchical_log_prob(params, config):
    lp = 0.0
    lp += np.sum(get_normal_lpdf(params['mu_k'],
                                 params['mu'],
                                 params['sd_mu']))
    lp += np.sum(get_normal_lpdf(params['tau_k'],
                                 params['tau'],
                                 params['sd_tau']))
    return lp


def get_observation_log_prob(params, data, config, weights=None):
    if weights is None:
        weights = np.ones(len(data['y']))
    lp = 0.0
    site_inds = data['site'] - 1 # site is 1-indexed, python 0-indexed
    loc = (params['mu_k'][site_inds] +
           params['tau_k'][site_inds] * data['treatment'])
    scale = params['sigma_k'][site_inds]

    use_smuber = data.get('use_smuber', False)
    if use_smuber:
        lp += np.sum(weights * get_smuber_lpdf(data['y'], loc, scale))
    else:
        lp += np.sum(weights * get_normal_lpdf(data['y'], loc, scale))

    return lp


def get_suff_observation_log_prob(params, data, config):
    lp = 0.0
    for k in range(data['num_sites']):
        beta = np.array([ params['mu_k'][k], params['tau_k'][k]])
        y_info = 1. / params['sigma_k'][k]**2
        lp += np.sum(get_regression_lpdf(
            xtx=data['xtx_k'][k],
            xty=data['xty_k'][k],
            yty=data['yty_k'][k],
            n_obs=data['n_obs_k'][k],
            beta=beta,
            y_info=y_info))

    return lp


def get_log_prob(params, data, config, prior_params, weights=None):
    use_prior = True
    if use_prior:
        log_prior = get_log_prior(params, prior_params, config)
    else:
        log_prior = 0.0
    hierarchical_log_prob = get_hierarchical_log_prob(params, config)
    observation_log_prob = get_observation_log_prob(params, data, config, weights)
    # observation_log_prob = get_suff_observation_log_prob(
    #     params, data, config)
    # print('\nDebugging')
    # print('{}\t{}\t{}\n'.format(
    #         log_prior,
    #         hierarchical_log_prob,
    #         observation_log_prob))
    return \
        log_prior + \
        hierarchical_log_prob + \
        observation_log_prob
