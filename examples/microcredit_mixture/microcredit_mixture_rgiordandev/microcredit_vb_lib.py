
import paragami
import vittles
import autograd
from autograd import numpy as np
from autograd import scipy as sp

def get_param_pattern(num_sites):
    num_sites = int(num_sites)
    param_pattern = paragami.PatternDict(free_default=True)

    # real mu[2];
    # real tau[2];
    # real<lower=0> sd_mu[2];
    # real<lower=0> sd_tau[2];
    param_pattern['mu'] = paragami.NumericVectorPattern(2)
    param_pattern['tau'] = paragami.NumericVectorPattern(2)
    param_pattern['sd_mu'] = paragami.NumericVectorPattern(2, lb=0.0)
    param_pattern['sd_tau'] = paragami.NumericVectorPattern(2, lb=0.0)

    # real sigma_control[2];
    # real sigma_TE[2];
    param_pattern['sigma_control'] = paragami.NumericVectorPattern(2)
    param_pattern['sigma_TE'] = paragami.NumericVectorPattern(2)

    # real<lower=0> sd_sigma_control[2];
    # real<lower=0> sd_sigma_TE[2];
    param_pattern['sd_sigma_control'] = paragami.NumericVectorPattern(2, lb=0.0)
    param_pattern['sd_sigma_TE'] = paragami.NumericVectorPattern(2, lb=0.0)

    # matrix[K,2] mu_k;
    # matrix[K,2] tau_k;
    # matrix[K,2] sigma_control_k;
    # matrix[K,2] sigma_TE_k;
    param_pattern['mu_k'] = paragami.NumericArrayPattern((num_sites, 2))
    param_pattern['tau_k'] = paragami.NumericArrayPattern((num_sites, 2))
    param_pattern['sigma_control_k'] = \
        paragami.NumericArrayPattern((num_sites, 2))
    param_pattern['sigma_TE_k'] = paragami.NumericArrayPattern((num_sites, 2))

    x_dim = 2 # P in Stan, used only for the categorical classifier
    num_categories = 3 # M in stan

    # matrix[M-1,P] beta; // the parent parameters minus the Mth category
    # matrix<lower=0>[M,P] sigma; // the set of M*P parent variances (not a covariance matrix)
    # matrix[M,P] beta_k_raw[K]; // the hierarchical increments
    param_pattern['beta'] = \
        paragami.NumericArrayPattern((num_categories - 1, x_dim))
    param_pattern['sigma'] = \
        paragami.NumericArrayPattern((num_categories, x_dim), lb=0.0)
    param_pattern['beta_k_raw'] = \
        paragami.NumericArrayPattern((num_sites, num_categories, x_dim))

    return param_pattern


def get_prior_param_pattern():
    prior_pattern = paragami.PatternDict()
    def append_loc_scale(param_name):
        prior_pattern[param_name + '_loc'] = paragami.NumericScalarPattern()
        prior_pattern[param_name + '_scale'] = \
            paragami.NumericScalarPattern(lb=0.0)

    append_loc_scale('beta')
    append_loc_scale('sigma')
    append_loc_scale('beta_k_raw')
    append_loc_scale('mu')
    append_loc_scale('tau')
    append_loc_scale('sd_mu')
    append_loc_scale('sd_tau')
    append_loc_scale('sigma_control')
    append_loc_scale('sd_sigma_control')
    append_loc_scale('sigma_TE')
    append_loc_scale('sd_sigma_TE')
    return prior_pattern


def get_default_prior_params():
    prior_params = dict()
    def append_loc_scale(param_name, loc, scale):
        prior_params[param_name + '_loc'] = loc
        prior_params[param_name + '_scale'] = scale

    append_loc_scale('beta', 0, 5)
    append_loc_scale('sigma', 0, 2)
    append_loc_scale('beta_k_raw', 0, 1)
    append_loc_scale('mu', 0, 100)
    append_loc_scale('tau', 0, 100)
    append_loc_scale('sd_mu', 0, 2)
    append_loc_scale('sd_tau', 0, 2)
    append_loc_scale('sigma_control', 0, 100)
    append_loc_scale('sd_sigma_control', 0, 2)
    append_loc_scale('sigma_TE', 0, 100)
    append_loc_scale('sd_sigma_TE', 0, 2)
    return prior_params


def get_normal_lpdf(z, mu, sigma):
    info = (1 / sigma**2)
    return -0.5 * info * (z - mu)**2 + 0.5 * np.log(info)

def get_cauchy_lpdf(z, loc, scale):
    z_norm = (z - loc) / scale
    return -np.log(scale) - np.log(1 + z_norm ** 2)

# Smuber = "Smooth Huber".  The -log(sigma) comes from the rescaling Jacobian.
def get_smuber_lpdf(z, mu, sigma):
    z_norm = (z - mu) / sigma
    z_norm2 = z_norm ** 2
    return -1 * z_norm2 / np.sqrt(1 + z_norm2) - np.log(sigma)

def get_log_prior(params, prior_params):
    # NOTE:
    # Parameters constrained to be positive but given unconstrained
    # symmetric priors need to be accounted for explicitly if the
    # means are allowed to be nonzero.

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
    # to_vector(beta) ~ normal(0,5);
    #lp += np.sum(get_normal_lpdf(params['beta'], 0, 5))
    lp += get_normal_prior('beta')

    # to_vector(sigma) ~ cauchy(0,2);
    #lp += np.sum(get_cauchy_lpdf(params['sigma'], 0, 2))
    lp += get_cauchy_prior('sigma')

    # for (m in 1:M){
    #   for (k in 1:K){
    #     beta_k_raw[k,m] ~ normal(0,1);
    #   }}
    #lp += np.sum(get_normal_lpdf(params['beta_k_raw'], 0, 1))
    lp += get_normal_prior('beta_k_raw')

    # mu ~ normal(0,100);
    # tau ~ normal(0,100);
    # lp += np.sum(get_normal_lpdf(params['mu'], 0, 100))
    # lp += np.sum(get_normal_lpdf(params['tau'], 0, 100))
    lp += get_normal_prior('mu')
    lp += get_normal_prior('tau')

    # sd_mu ~ cauchy(0,2);
    # sd_tau ~ cauchy(0,2);
    # lp += np.sum(get_cauchy_lpdf(params['sd_mu'], 0, 2))
    # lp += np.sum(get_cauchy_lpdf(params['sd_tau'], 0, 2))
    lp += get_cauchy_prior('sd_mu')
    lp += get_cauchy_prior('sd_tau')

    # sigma_control ~ normal(0,100);
    #lp += np.sum(get_normal_lpdf(params['sigma_control'], 0, 100))
    lp += get_normal_prior('sigma_control')

    # sd_sigma_control ~ cauchy(0,2);
    #lp += np.sum(get_cauchy_lpdf(params['sd_sigma_control'], 0, 2))
    lp += get_cauchy_prior('sd_sigma_control')

    # sigma_TE ~ normal(0,100);
    #lp += np.sum(get_normal_lpdf(params['sigma_TE'], 0, 100))
    lp += get_normal_prior('sigma_TE')

    # sd_sigma_TE ~ cauchy(0,2);
    #lp += np.sum(get_cauchy_lpdf(params['sd_sigma_TE'], 0, 2))
    lp += get_cauchy_prior('sd_sigma_TE')

    return lp


def get_hierarchical_log_prob(params):
    # for (k in 1:K){ // Site
    #   for (i in 1:2){ // Sign
    #     mu_k[k,i] ~ normal(mu[i], sd_mu[i]);
    #     tau_k[k,i] ~ normal(tau[i], sd_tau[i]);
    #     sigma_control_k[k,i] ~ normal(sigma_control[i], sd_sigma_control[i]);
    #     sigma_TE_k[k,i] ~ normal(sigma_TE[i], sd_sigma_TE[i]);
    #   }
    # }

    lp = 0.0
    lp += np.sum(get_normal_lpdf(params['mu_k'],
                                 params['mu'][None, :],
                                 params['sd_mu'][None, :]))
    lp += np.sum(get_normal_lpdf(params['tau_k'],
                                 params['tau'][None, :],
                                 params['sd_tau'][None, :]))
    lp += np.sum(get_normal_lpdf(params['sigma_control_k'],
                                 params['sigma_control'][None, :],
                                 params['sd_sigma_control'][None, :]))
    lp += np.sum(get_normal_lpdf(params['sigma_TE_k'],
                                 params['sigma_TE'][None, :],
                                 params['sd_sigma_TE'][None, :]))
    return lp


def get_log_softmax(log_odds):
    # See math/stan/math/prim/mat/prob/categorical_logit_lpmf.hpp
    # and stan/math/prim/mat/fun/log_softmax.hpp, though due to the special
    # structure of this problem I'll do categorical_logit_lpmf by hand.
    log_probs = log_odds - sp.special.logsumexp(log_odds, axis=1)[:, None]
    return log_probs


def get_weights(weights, data, normalize=False):
    num_obs = len(data['profit'])
    if weights is None:
        w = np.ones(num_obs)
    else:
        w = weights
    if normalize:
        w = w / np.sum(w)
    return w


def get_category_log_prob(params, data, weights=None):
    lp = 0.0
    weights = get_weights(weights, data)

    # beta_full = append_row(beta, rep_row_vector(0, P));
    x_dim = params['beta'].shape[1]
    beta_full = np.vstack([params['beta'], np.zeros((1, x_dim))])
    #beta_full = np.hstack([params['beta'], np.zeros((x_dim, 1))]).T

    # for (m in 1:M){
    #  for (k in 1:K){
    #   for (p in 1:P){
    #    beta_k[k,m,p] = beta_full[m,p] + sigma[m,p] * beta_k_raw[k,m,p];
    # }}}
    beta_k = beta_full[None, :, :] + \
             params['sigma'][None, :, :] * params['beta_k_raw']

    # > X <- cbind(rep(1,N), data$treatment)
    # for (n in 1:N)
    #   cat[n] ~ categorical_logit(beta_k[site[n]] * x[n]);

    # The Stan implementation seems computationally wasteful, since we know
    # a priori that there are only two possible values for x.
    assert np.all(np.unique(data['treatment']) == np.array([0, 1]))
    assert np.all(np.unique(data['cat']) == np.array([1, 2, 3]))

    unique_sites = np.unique(data['site']).astype(int)
    unique_cats = np.unique(data['cat']).astype(int)

    #print(beta_k.shape) # Site, outcome, x.
    for treat in [0, 1]:
        data_rows = data['treatment'].astype(int) == treat
        x = np.array([1, treat])

        # These are arrays of (site * outcome).
        log_odds = np.einsum('kmp,p->km', beta_k, x)
        log_probs = get_log_softmax(log_odds)

        # If this is too slow, you could use grouped_sum.
        for site_i in range(len(unique_sites)):
            for cat_i in range(len(unique_cats)):
                site = int(unique_sites[site_i])
                cat = int(unique_cats[cat_i])
                row_count = np.sum(weights[data_rows] * \
                    np.logical_and(
                       data['site'][data_rows] == site,
                       data['cat'][data_rows] == cat))
                lp += row_count * log_probs[site - 1, cat - 1]

    return lp


def get_observation_log_prob(params, data, weights=None):
    # for (n in 1:N_neg){
    #     y_neg[n] ~ lognormal(
    #       mu_k[site_neg[n],1] + tau_k[site_neg[n],1] * treatment_neg[n],
    #       exp(sigma_control_k[site_neg[n],1] +
    #           sigma_TE_k[site_neg[n],1] * treatment_neg[n]));
    # }
    # for(n in 1:N_pos){
    #     y_pos[n] ~ lognormal(
    #       mu_k[site_pos[n],2] + tau_k[site_pos[n],2] * treatment_pos[n],
    #       exp(sigma_control_k[site_pos[n],2] +
    #           sigma_TE_k[site_pos[n],2] * treatment_pos[n]));
    # }

    lp = 0.0
    cat = data['cat']
    weights = get_weights(weights, data)

    # The entries of obs_categories corresponds to the columns of mu_k, &c.
    obs_categories = [1, 3]
    for i in range(len(obs_categories)):
        rows = (cat == obs_categories[i])
        site = data['site'][rows] - 1
        treat = data['treatment'][rows]

        loc = params['mu_k'][site, i] + params['tau_k'][site, i] * treat
        scale = np.exp(params['sigma_control_k'][site, i] +
                       params['sigma_TE_k'][site, i] * treat)

        log_abs_y = np.log(np.abs(data['profit'][rows]))
        use_smuber = data.get('use_smuber', False)
        if use_smuber:
            lp += np.sum(weights[rows] * get_smuber_lpdf(log_abs_y, loc, scale))
        else:
            lp += np.sum(weights[rows] * get_normal_lpdf(log_abs_y, loc, scale))

    return lp


def get_log_prob(params, data, prior_params, weights=None):
    log_prior = get_log_prior(params, prior_params)
    hierarchical_log_prob = get_hierarchical_log_prob(params)
    category_log_prob = get_category_log_prob(params, data, weights)
    observation_log_prob = get_observation_log_prob(params, data, weights)
    # print('\nDebugging')
    # print('{}\t{}\t{}\t{}\n'.format(
    #         log_prior,
    #         hierarchical_log_prob,
    #         category_log_prob,
    #         observation_log_prob))
    return \
        log_prior + \
        hierarchical_log_prob + \
        category_log_prob + \
        observation_log_prob


def get_log_prob_with_jacobian(params, data, prior_params, param_pattern,
                               weights=None):
    # get_log_prob returns a density in the folded space, so we need to
    # add log(det(dfolded / dfree))
    # jac = param_pattern.unfreeing_jacobian(params)
    log_det_jac = param_pattern.log_abs_det_unfreeing_jacobian(params)
    return get_log_prob(params, data, prior_params, weights) + log_det_jac


#################################
# ADVI


# TODO: use advi_lib

def get_advi_pattern(param_pattern):
    param_dim = int(param_pattern.flat_length(free=True))
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


def get_advi_loss(advi_params, base_draws, data, prior_params, param_pattern,
                  weights=None):
    param_free_mat = encode_draws(advi_params, base_draws)
    num_draws = base_draws.shape[1]
    lp_sum = 0.0
    for d in range(num_draws):
        params = param_pattern.fold(param_free_mat[:, d], free=True)
        lp_sum += get_log_prob_with_jacobian(
            params, data, prior_params, param_pattern, weights)

    advi_entropy = get_advi_entropy(advi_params)
    return -1 * (lp_sum / num_draws + advi_entropy)


def get_base_draws(num_draws, param_pattern):
    param_dim = param_pattern.flat_length(free=True)
    return np.random.normal(size=(param_dim, num_draws))


def get_advi_objective(base_draws, data, prior_params, advi_param_pattern,
                       param_pattern, weights=None):
    get_log_loss_free = paragami.FlattenFunctionInput(
        lambda advi_params: \
            get_advi_loss(advi_params, base_draws, data, prior_params,
                          param_pattern, weights),
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
