# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:55:46 2022

@author: cosmi
"""

import numpy as np
import cv2

from qutip import Qobj, rand_dm, fidelity, displace, qdiags, qeye, expect
from qutip.states import coherent, coherent_dm, thermal_dm, fock_dm
from qutip.visualization import plot_wigner, hinton, matrix_histogram
from qutip.wigner import qfunc
import qutip

import matplotlib.pyplot as plt
from matplotlib import animation

# some pretty printing and animation stuff
from IPython.display import clear_output

"""
Iterative Maximum Likelihood estimation based on photon number counting.
"""


def measure_population(beta, rho):
    """
    Measures the photon number statistics for state rho when displaced
    by angle alpha.
    
    Parameters
    ----------    
    alpha: np.complex
        A complex displacement.

    rho:
        The density matrix as a QuTiP Qobj (`qutip.Qobj`)

    Returns
    -------
    population: ndarray
        A 1D array for the probabilities for populations.
    """
    hilbertsize = rho.shape[0]
    # Apply a displacement to the state and then measure the diagonals.

    D = displace(hilbertsize, beta)
    rho_disp = D*rho*D.dag()
    populations = np.real(np.diagonal(rho_disp.full()))
    return populations

def roperator(beta, rho, measured):
    """
    Calculates the iterative ratio operator for measured probability for photons
    (j) to the analytical prediction for some rho.

    Parameters
    ----------
    beta: list_like
        A list of the displacements that were applied to the state before
        measurement.

    rho: `qutip.Qobj`
        The current estimate of the density matrix.

    measured: list_like
        A list of the measurement statistics (diagonal terms) for each beta.

    Returns
    -------
    R: `qutip.Qobj`
        The iterative operator which we are going to apply for state
        reconstruction.
    """
    hilbert_size = rho.shape[0]

    # initialize an empty operator and build it

    R = 0*qeye(hilbert_size)
    calculated_measurements = measure_population(beta, rho)

    for n in range(hilbert_size):
        op = fock_dm(hilbert_size, n)
        D = displace(hilbert_size, beta)
        displaced_D = D.dag()*op*D
        ratio = measured[n]/(calculated_measurements[n] + 1e-6)
        displaced_D = ratio*displaced_D
        R += displaced_D

    return R


hilbert_size = 4
alpha_range = 1.9

alphas = np.array([alpha_range, -alpha_range - 1j*alpha_range,
                   -alpha_range + 1j*alpha_range])

rho_true = sum([coherent_dm(hilbert_size, a) for a in alphas])/3


betas = [1.7, -2, 2.2j, -2.1 - 2.4j, -2 + 2j]
measured_populations = [measure_population(b, rho_true) for b in betas]
width = 1

random_rho = rand_dm(hilbert_size)
hinton(random_rho)
plt.show()



# visualize H

lbls_list = [[str(d) for d in range(hilbert_size)], ["u", "d"]]

xlabels = []

fig, ax = matrix_histogram(random_rho, xlabels, xlabels, limits=[-1,1])

ax.view_init(azim=-55, elev=45)

plt.show()




fig, ax = plt.subplots(1, 3, figsize=(15, 5))
indices = np.arange(hilbert_size)

plot_wigner(random_rho, fig, ax[0])
ax[0].scatter(np.real(betas), np.imag(betas), marker="x")
ax[0].set_title("Random inital state wigner function")
for i in range(len(betas)):
    ax[1].bar(indices, measured_populations[i],
              label = r"$beta = {}$".format(i), width=(i+1)/12)

ax[1].set_title("Population measurement statistics")
ax[1].set_xlabel("n")
ax[1].set_ylabel("Photon number probability")


plot_wigner(rho_true, fig, ax[2])
ax[2].scatter(np.real(betas), np.imag(betas), marker="x")
ax[2].set_title("Target state wigner function")


plt.savefig("Target state wigner function.png")

plt.show()
