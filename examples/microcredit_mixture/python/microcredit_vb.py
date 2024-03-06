#!/usr/bin/env python
# Example invocation:
# ./microcredit_vb.py --git_repo_loc=/home/rgiordan/Documents/git_repos/AdversarialInfluenceWorkbench --num_draws=20

# Before running, follow the installation instructions in the `AdversarialInfluenceWorkbench/examples/microcredit_mixture/README.md` file.

import argparse
import paragami
#import vittles
#import autograd
import numpy as np
#from autograd import numpy as np
#from autograd import scipy as sp
#import matplotlib.pyplot as plt
import os
import pandas as pd
import time
import scipy as osp


import microcredit_mixture_rgiordandev
from microcredit_mixture_rgiordandev import microcredit_vb_lib as mm_lib
from microcredit_mixture_rgiordandev import \
    get_log_prob, get_param_pattern, get_log_prob_with_jacobian

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

np.random.seed(args.seed)
num_draws = args.num_draws

# This controls whether we use a "smoothed Huber" loss rather
# than a squared loss to try to mitigate outlier sensitivity.  
# It was experimental and we don't use it in the paper.
use_smuber = False


raw_data = pd.read_csv(
    os.path.join(args.git_repo_loc,
                 'examples/microcredit_mixture/R/data',
                 'microcredit_project_data_cleaned.csv'))

data = dict()

data['site'] = raw_data["site"].to_numpy().astype('int')
data['profit'] = raw_data["profit"].to_numpy()
data['treatment'] = raw_data["treatment"].to_numpy().astype('int')
data['cat'] = raw_data["cat"].to_numpy().astype('int')
data['use_smuber'] = use_smuber

prior_param_pattern = mm_lib.get_prior_param_pattern()
prior_params = mm_lib.get_default_prior_params()
prior_param_pattern.flatten(prior_params, free=False)


assert np.min(data['site']) == 1
num_sites = int(np.max(data['site'])) # K in stan.  Convert to int to allow json serialization.

param_pattern = get_param_pattern(num_sites)
params_rand = param_pattern.random()


advi_param_pattern = mm_lib.get_advi_pattern(param_pattern)
base_draws = mm_lib.get_base_draws(num_draws, param_pattern)

advi_params = advi_param_pattern.random()
#mm_lib.get_advi_loss(advi_params, base_draws, data, prior_params, param_pattern)



# Get the DADVI draws and define an objective

base_draws = mm_lib.get_base_draws(num_draws, param_pattern)
objective = mm_lib.get_advi_objective(
    base_draws, data, prior_params, advi_param_pattern, param_pattern)

advi_params_rand = advi_param_pattern.random()
advi_params_rand_free = advi_param_pattern.flatten(advi_params_rand)


# ## Run an initial fit.

print('Running initial optimization.')

init_advi_params = advi_param_pattern.empty(valid=True)
init_advi_params['mean'][:] = 0.0
init_advi_params['log_sd'][:] = -1

opt_init_time = time.time()
objective.reset()
objective.set_log_every(1)
opt_result = osp.optimize.minimize(
    method='trust-ncg',
    x0=advi_param_pattern.flatten(init_advi_params, free=True),
    fun=objective.f,
    jac=objective.grad,
    hessp=objective.hessian_vector_product,
    options={'maxiter': 100})
opt_init_time = time.time() - opt_init_time
print(f'Initial optimization took f{opt_init_time / 60} minutes')


# ## Preconditioned optimization

# Evaluate the Hessian at the initial fit, check for positive definiteness, and then re-run with preconditioning to make sure we're at a good optimum.

print('Running preconditioned optimization.')

print('Computing Hessian.')
hess_pc = objective.hessian(opt_result.x)
evs = np.linalg.eigvals(hess_pc)
#np.min(evs), np.max(evs)

get_log_loss_free = paragami.FlattenFunctionInput(
    lambda advi_params: \
        mm_lib.get_advi_loss(advi_params, base_draws, data,
                             prior_params, param_pattern, weights=None),
    patterns=advi_param_pattern,
    free=True)

advi_loss_precond = paragami.PreconditionedFunction(get_log_loss_free)
advi_loss_precond.set_preconditioner_with_hessian(hessian=hess_pc)
objective_precond = paragami.OptimizationObjective(advi_loss_precond)
init_x = advi_loss_precond.precondition(opt_result.x)

print('Optimizing.')
opt_pc_time = time.time()
objective_precond.reset()
objective_precond.set_log_every(1)
opt_result = osp.optimize.minimize(
    method='trust-ncg',
    x0=init_x,
    fun=objective_precond.f,
    jac=objective_precond.grad,
    hessp=objective_precond.hessian_vector_product,
    options={'maxiter': 500})
print(opt_result.message)
opt_pc_time = time.time() - opt_pc_time
print(f'Preconditioned optimization took f{opt_pc_time / 60} minutes')


# The total optimization time is the sum of both.
opt_time = opt_init_time + opt_pc_time

print(f'Total optimization took f{opt_time / 60} minutes')


# ## Inspect and save results.

advi_params = advi_param_pattern.fold(advi_loss_precond.unprecondition(opt_result.x))
advi_params_free = advi_param_pattern.flatten(advi_params, free=True)
params_opt = param_pattern.fold(advi_params['mean'])

hess_time = time.time()
hess0 = objective.hessian(advi_params_free)
hess_time = time.time() - hess_time
evs = np.linalg.eigvals(hess0)

print(f'Hessian computation time: f{hess_time / 60} minutes')
print(f'Max ev: {np.max(evs)}, Min ev: {np.min(evs)}, Condition number: {np.max(evs) / np.min(evs)}')


grad0 = objective.grad(advi_params_free)

print(f'Gradient norm: {np.linalg.norm(grad0)}')
print(f'Newton step norm: {np.linalg.norm(np.linalg.solve(hess0, grad0))}')


# In[56]:


save_dict = {}
save_dict['opt_time'] = opt_time
save_dict['hess_time'] = hess_time
save_dict['hess0'] = hess0
save_dict['advi_params_free'] = advi_param_pattern.flatten(advi_params, free=True)
save_dict['params_free'] = param_pattern.flatten(params_opt, free=True)
save_dict['advi_param_pattern_json'] = advi_param_pattern.to_json()
save_dict['param_pattern_json'] = param_pattern.to_json()
save_dict['base_draws'] = base_draws
save_dict['use_smuber'] = data['use_smuber']
#save_dict['optimization_log'] = objective_precond.optimization_log # No longer works with npz
save_dict['prior_param_pattern_json'] = prior_param_pattern.to_json()
save_dict['prior_params_flat'] = prior_param_pattern.flatten(prior_params, free=False)


# In[57]:
out_filename = 'microcredit_project_advi_{}draws_smuber{}.npz'.format(
    args.num_draws, data['use_smuber'])

print(f'Saving to {os.path.join(out_path, out_filename)}')
np.savez_compressed(**save_dict, file=out_filename)

