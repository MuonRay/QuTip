# -*- coding: utf-8 -*-
"""
Created on Fri Dec 23 02:16:32 2022

@author: cosmi
"""

import numpy as np
import qutip as qt
import qinfer as qi
import matplotlib
matplotlib.style.use('ggplot')
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse


from qutip.ipynbtools import plot_animation

from mpl_toolkits.axes_grid1 import make_axes_locatable

import tempfile as tf
import os

plt.rcParams['legend.fontsize'] = 18
plt.rcParams['axes.labelsize'] = 22
from qinfer.tomography.plotting_tools import (
    plot_cov_ellipse, plot_decorate_rebits, plot_rebit_modelparams, plot_rebit_prior, plot_rebit_posterior
)
I, X, Y, Z = qt.qeye(2), qt.sigmax(), qt.sigmay(), qt.sigmaz()


n_steps=1000 
save_stills=()
save_stills=[1000 * idx for idx in range(1, 10)]

tempdir = tf.mkdtemp()
print("Saving to {}.".format(tempdir))
pbar = qt.ui.TextProgressBar()
pbar.start(n_steps, chunk_size=1)

basis = qi.tomography.pauli_basis(1)
diffusive_model = qi.BinomialModel(qi.tomography.DiffusiveTomographyModel(basis, False))
    
true_states = []
est_states = []
    
true_state = np.array([[1, 0.99, 0, 0, 0]]) / np.sqrt(2)
true_state[0, -1] = 0.0045
diffusive_prior = qi.ProductDistribution(
qi.tomography.GADFLIDistribution(
qi.tomography.GinibreReditDistribution(basis),
I / 2 - 0.02 * Z / 2 + 0.88 * X / 2),
qi.LogNormalDistribution(0.0, 0.006)
)
    
updater = qi.smc.SMCUpdater(diffusive_model, 3000, diffusive_prior)
heuristic = qi.tomography.RandomPauliHeuristic(updater, other_fields={'n_meas': 25, 't': 1})
    
fig, (ax_rebit, ax_cov) = plt.subplots(1, 2, figsize=(1920 * 12 / 1080, 12))




for idx_exp in range(n_steps):
    expparams = heuristic()
    outcome = diffusive_model.simulate_experiment(true_state, expparams)
    updater.update(outcome, expparams)
    est_states.append(basis.modelparams_to_state(updater.est_mean()[:-1]))
    true_states.append(basis.modelparams_to_state(true_state[0, :-1]))
    err = np.linalg.norm(true_state[0] - updater.est_mean(), 2)
        
    true_state = diffusive_model.update_timestep(true_state, expparams)[:, :, 0]
        
    plt.sca(ax_rebit)
    plot_rebit_posterior(updater, prior=None, true_state=true_state, rebit_axes=[1, 3], true_size=600)
    plt.xticks([-1, 0, 1], size=14)
    plt.yticks([-1, 0, 1], size=14)
        
    plt.sca(ax_cov)
    updater.plot_covariance(param_slice=np.s_[1:4], tick_params={'size': 20})
        
    plt.savefig(os.path.join(tempdir, '{:05}.png'.format(idx_exp)), dpi=1080 * 100 / 1200)
    
    #plot_animation(fig, outcome, true_state)


#ani = animation.FuncAnimation(fig, diffusion_movie, frames=30)
#FFwriter = animation.FFMpegWriter()
#ani.save('plot.mp4', writer=FFwriter)
