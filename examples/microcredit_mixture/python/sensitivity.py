#!/usr/bin/env python
# Example invocation:
# ./sensitivity.py --git_repo_loc=/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench --num_draws=30


import paragami
import vittles
import autograd
from autograd import numpy as np
from autograd import scipy as sp
#import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
from copy import deepcopy
import time
import scipy as osp

import microcredit_mixture_rgiordandev
from microcredit_mixture_rgiordandev import microcredit_vb_lib as mm_lib

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('--seed', type=int, help='RNG seed', default=42)
parser.add_argument('--num_draws', type=int, help='RNG seed', default=30)
parser.add_argument('--git_repo_loc', type=str, help='Path to git repository', default=None)

args = parser.parse_args()

if not args.git_repo_loc:
    raise ValueError('You must specify --git_repo_loc')

if args.num_draws <= 1:
    raise ValueError('num_draws must be > 1')

out_path = os.path.join(
    args.git_repo_loc,
    'examples/microcredit_mixture/python/output')
if not os.path.exists(out_path):
    raise ValueError(f'Output path {out_path} does not exist')


num_draws = args.num_draws
use_smuber = False


# Load the initial fit and data
out_path = os.path.join(
    args.git_repo_loc,
    'examples/microcredit_mixture/python/output')

initial_fit_filename = \
    os.path.join(out_path, f'microcredit_project_advi_{num_draws}draws_smuber{use_smuber}.npz')

#out_filename = 'output/microcredit_project_advi_{}draws.npz'.format(5)
advi_dict = np.load(initial_fit_filename)

advi_param_pattern = paragami.get_pattern_from_json(str(advi_dict['advi_param_pattern_json']))
param_pattern = paragami.get_pattern_from_json(str(advi_dict['param_pattern_json']))
advi_param_pattern = paragami.get_pattern_from_json(str(advi_dict['advi_param_pattern_json']))
param_pattern = paragami.get_pattern_from_json(str(advi_dict['param_pattern_json']))
hess0 = advi_dict['hess0']
advi_params_free = advi_dict['advi_params_free']
advi_params = advi_param_pattern.fold(advi_params_free, free=True)
mean_params = param_pattern.fold(advi_dict['params_free'], free=True)
base_draws = advi_dict['base_draws']
prior_pattern = paragami.get_pattern_from_json(str(advi_dict['prior_param_pattern_json']))
prior_params = prior_pattern.fold(advi_dict['prior_params_flat'], free=False)


# Load the raw data

raw_data = pd.read_csv(
    os.path.join(
        args.git_repo_loc,
        'examples/microcredit_mixture/R/data/microcredit_project_data_cleaned.csv'))

data = dict()
data['site'] = raw_data["site"].to_numpy().astype('int')
data['profit'] = raw_data["profit"].to_numpy()
data['treatment'] = raw_data["treatment"].to_numpy().astype('int')
data['cat'] = raw_data["cat"].to_numpy().astype('int')
data['weights'] = None

num_obs = len(data['profit'])
data['use_smuber'] = use_smuber

# Define the loss function

get_log_loss_free = paragami.FlattenFunctionInput(
    lambda advi_params: mm_lib.get_advi_loss(
        advi_params, base_draws, data, prior_params, param_pattern),
    patterns=advi_param_pattern,
    free=True)
objective = paragami.OptimizationObjective(get_log_loss_free)
grad0 = objective.grad(advi_params_free)

print('Gradient and netwon step norm at optimum (should be close to zero):')
print(np.linalg.norm(grad0))
print(np.linalg.norm(np.linalg.solve(hess0, grad0)))


def get_weighted_log_loss(advi_params, weights):
    return mm_lib.get_advi_loss(
        advi_params, base_draws, data, prior_params, param_pattern, weights)

get_weighted_log_loss_free = paragami.FlattenFunctionInput(
    get_weighted_log_loss,
    patterns=advi_param_pattern,
    free=True)

estimating_equation = autograd.grad(get_weighted_log_loss, argnum=0)
w1 = np.ones(num_obs)

# Sanity checks
assert     get_weighted_log_loss_free(advi_params_free, w1) -         get_log_loss_free(advi_params_free) == 0

assert     get_weighted_log_loss_free(advi_params_free, None) -         get_log_loss_free(advi_params_free) == 0


# Make sure this is feasible.
w_grad = autograd.grad(get_weighted_log_loss_free, argnum=1)(advi_params_free, w1)
print(w_grad.shape)


get_theta_grad = autograd.grad(get_weighted_log_loss_free, argnum=0)
get_cross_hess = autograd.jacobian(get_theta_grad, argnum=1)


print('Computing cross Hessian:')
cross_hess_time = time.time()
cross_hess_at_opt = get_cross_hess(advi_params_free, w1)
cross_hess_time = time.time() - cross_hess_time


print(cross_hess_time)
print(cross_hess_at_opt.shape)


print('Computing sensitivity.')
hess_solver = vittles.solver_lib.get_dense_cholesky_solver(hess0)

weight_sens_approx =     vittles.sensitivity_lib.HyperparameterSensitivityLinearApproximation(
        objective_fun=get_weighted_log_loss_free,
        opt_par_value=advi_params_free,
        hyper_par_value=w1,
        validate_optimum=False,
        hessian_at_opt=hess0,
        cross_hess_at_opt=cross_hess_at_opt)


influence_mat = weight_sens_approx.get_dopt_dhyper()

print(influence_mat.shape)


print('Computing LR covariance')
#get_advi_mean(advi_params, base_draws, param_pattern)
get_advi_mean_free = paragami.FlattenFunctionInput(
    mm_lib.get_advi_mean, advi_param_pattern, free=True, argnums=[0])
advi_mean_free_jacobian = autograd.jacobian(
    get_advi_mean_free, argnum=0)(advi_params_free, base_draws, param_pattern)
mean_influence_mat = advi_mean_free_jacobian @ influence_mat


lr_cov = advi_mean_free_jacobian @ hess_solver(advi_mean_free_jacobian.T)

def get_advi_sd(advi_params, base_draws, param_pattern):
    flat_draws = mm_lib.get_flat_draws(advi_params, base_draws, param_pattern)
    return np.std(flat_draws, axis=0)

advi_sd = get_advi_sd(advi_params, base_draws, param_pattern)
lrvb_sd = np.sqrt(np.diag(lr_cov))
advi_mean_free_jacobian.shape


# Save results

save_dict = {}
save_dict['out_filename'] = initial_fit_filename
save_dict['cross_hess_time'] = cross_hess_time
save_dict['cross_hess_at_opt'] = cross_hess_at_opt
save_dict['advi_mean_free_jacobian'] = advi_mean_free_jacobian
save_dict['influence_mat'] = influence_mat
save_dict['lr_cov'] = lr_cov
save_dict['mean_influence_mat'] = mean_influence_mat
save_dict['use_smuber'] = data['use_smuber']


sens_filename = os.path.join(
    out_path,
    f'microcredit_project_weight_sensitivity_{num_draws}draws_smuber{use_smuber}.npz')
print(f'Saving to {sens_filename}')

np.savez_compressed(**save_dict, file=sens_filename)
print('Success!')