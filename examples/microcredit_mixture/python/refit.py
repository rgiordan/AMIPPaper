#!/usr/bin/env python

"""Load a previous fit and refit, possibly with new weights.

This uses the output of the jupyter notebook microcredit_vb.

Example usage:

./refit.py \
    --initial_fit_filename=
    --data_filename=
    --reweight_filename=
    --output_filename=
"""

import argparse

import paragami
import vittles
import autograd
from autograd import numpy as np
import pandas as pd
import scipy as sp
import os

import time

import microcredit_mixture_rgiordandev
from microcredit_mixture_rgiordandev import microcredit_vb_lib as mm_lib

parser = argparse.ArgumentParser()

parser.add_argument('--seed', default=42, type=int)
parser.add_argument('--initial_fit_filename', default=None, type=str)
parser.add_argument('--data_filename', default=None, type=str)
parser.add_argument('--reweight_filename', default=None, type=str)
parser.add_argument('--output_filename', default=None, type=str)

args = parser.parse_args()

np.random.seed(args.seed)

output_dir, _  = os.path.split(args.output_filename)
if not os.path.isdir(output_dir):
    raise ValueError('Output directory {} does not exist'.format(output_dir))

if not os.path.isfile(args.initial_fit_filename):
    raise ValueError('Initial fit file {} does not exist'.format(
        args.initial_fit_filename))

if not os.path.isfile(args.data_filename):
    raise ValueError('Data file {} does not exist'.format(
        args.data_filename))

if not os.path.isfile(args.reweight_filename):
    raise ValueError('Weight file {} does not exist'.format(
        args.reweight_filename))

##################
# Load and run


# data
raw_data = pd.read_csv(args.data_filename)
data = dict()
data['site'] = raw_data["site"].to_numpy().astype('int')
data['profit'] = raw_data["profit"].to_numpy()
data['treatment'] = raw_data["treatment"].to_numpy().astype('int')
data['cat'] = raw_data["cat"].to_numpy().astype('int')
n_obs = len(data['profit'])

# Load initial fit
advi_dict = np.load(args.initial_fit_filename)
advi_param_pattern = paragami.get_pattern_from_json(
    str(advi_dict['advi_param_pattern_json']))
param_pattern = paragami.get_pattern_from_json(
    str(advi_dict['param_pattern_json']))
hess0 = advi_dict['hess0']
advi_params_free = advi_dict['advi_params_free']
advi_params = advi_param_pattern.fold(advi_params_free, free=True)
base_draws = advi_dict['base_draws']
prior_param_pattern = paragami.get_pattern_from_json(
    str(advi_dict['prior_param_pattern_json']))
prior_params = prior_param_pattern.fold(
    advi_dict['prior_params_flat'], free=False)

reweight_dict = np.load(args.reweight_filename)
weights = reweight_dict['weights']
print('Leaving out {} observations'.format(np.sum(1 - weights)))

# Define objective
get_log_loss_free = paragami.FlattenFunctionInput(
    lambda advi_params: \
        mm_lib.get_advi_loss(advi_params, base_draws, data,
                             prior_params, param_pattern, weights),
    patterns=advi_param_pattern,
    free=True)

advi_loss_precond = paragami.PreconditionedFunction(get_log_loss_free)
advi_loss_precond.set_preconditioner_with_hessian(hessian=hess0)
objective_precond = paragami.OptimizationObjective(advi_loss_precond)
init_x = advi_loss_precond.precondition(advi_params_free)

# Optimize
print('Optimizing!  -ㅅ-')
opt_time = time.time()
objective_precond.reset()
objective_precond.set_log_every(1)
opt_result = sp.optimize.minimize(
    method='trust-ncg',
    x0=init_x,
    fun=objective_precond.f,
    jac=objective_precond.grad,
    hessp=objective_precond.hessian_vector_product,
    options={'maxiter': 500})
opt_time = time.time() - opt_time

print(opt_result.message)
print('Optimized!  ^-^')
print('Optimization took {} seconds'.format(opt_time))

advi_params_refit = advi_param_pattern.fold(
    advi_loss_precond.unprecondition(opt_result.x), free=True)
advi_params_refit_free = advi_param_pattern.flatten(advi_params_refit, free=True)
params_opt_refit = param_pattern.fold(advi_params['mean'], free=True)

# Get the LR covariance
print('Computing Hessian')
objective = mm_lib.get_advi_objective(
    base_draws, data, prior_params, advi_param_pattern, param_pattern,
    weights=weights)

hess_time = time.time()
hess1 = objective.hessian(advi_params_refit_free)
hess_time = time.time() - hess_time
evs = np.linalg.eigvals(hess1)

print(hess_time)
print(np.max(evs), np.min(evs), np.max(evs) / np.min(evs))

print('Computing LR covariance')
hess_solver = vittles.solver_lib.get_dense_cholesky_solver(hess1)

get_advi_mean_free = paragami.FlattenFunctionInput(
    mm_lib.get_advi_mean, advi_param_pattern, free=True, argnums=[0])
advi_mean_free_jacobian = autograd.jacobian(
    get_advi_mean_free, argnum=0)(advi_params_refit_free, base_draws, param_pattern)
lr_cov = advi_mean_free_jacobian @ hess_solver(advi_mean_free_jacobian.T)

# Save the results

save_dict = {}
save_dict['opt_time'] = opt_time
save_dict['initial_fit_filename'] = args.initial_fit_filename
save_dict['reweight_filename'] = args.reweight_filename
save_dict['weights'] = weights
save_dict['lr_cov'] = lr_cov
save_dict['advi_params_free'] = \
    advi_param_pattern.flatten(advi_params_refit, free=True)
save_dict['advi_param_pattern_json'] = advi_param_pattern.to_json()
save_dict['param_pattern_json'] = param_pattern.to_json()
save_dict['base_draws'] = base_draws

print('Saving to {}'.format(args.output_filename))
np.savez_compressed(**save_dict, file=args.output_filename)
print('Done!  ◕‿◕')
