# -*- coding: utf-8 -*-
"""
Created on Fri May  5 14:35:27 2023

@author: cosmi
"""

# imports 

import numpy as np

from qutip import Qobj, rand_dm, fidelity, displace, qdiags, qeye, expect
from qutip.states import coherent, coherent_dm, thermal_dm, fock_dm
from qutip.random_objects import rand_dm
from qutip.visualization import plot_wigner, hinton, plot_wigner_fock_distribution
from qutip.wigner import qfunc
import qutip

import matplotlib.pyplot as plt
from matplotlib import animation, colors


from IPython.display import clear_output

hilbert_size = 32
psi = coherent(hilbert_size, 0)
d = displace(hilbert_size, 2+2j)

fig, ax = plt.subplots(1, 4, figsize=(19, 4))

plot_wigner_fock_distribution(psi, fig=fig, axes=[ax[0], ax[1]])
plot_wigner_fock_distribution(d*psi, fig=fig, axes=[ax[2], ax[3]])

ax[0].set_title(r"Initial state, $\psi_{vac} = |0 \rangle$")
ax[2].set_title(r"Displaced state, $D(\alpha=2+2i )\psi_{vac}$")
plt.show()

alpha_range = 2
alphas = np.array([alpha_range, -alpha_range - 1j*alpha_range,
                   -alpha_range + 1j*alpha_range])

psi = sum([coherent(hilbert_size, a) for a in alphas])
psi = psi.unit()
rho = psi*psi.dag()


fig, ax = plot_wigner_fock_distribution(rho, figsize=(9, 4))
ax[0].set_title("Superposition of three coherent states")
plt.show()


def measure_q(beta, rho):
    """
    Measures the generalized q function values for the state density matrix.
    
    Parameters
    ----------    
    beta: np.complex
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

    D = displace(hilbertsize, -beta)
    rho_disp = D*rho*D.dag()
    
    # measure all the elements in the diagonal
    populations = np.real(np.diagonal(rho_disp.full()))
    return populations


betas = [1.7, -2, 2.5j, -2.1 - 2.1j, -2 + 2j]
generalized_Q = [measure_q(b, rho) for b in betas]


fig, ax = plt.subplots(1, 3, figsize=(15, 4))
indices = np.arange(hilbert_size)

plot_wigner(rho, fig, ax[0])
ax[0].scatter(np.real(betas), np.imag(betas), marker="x")
ax[0].set_title(r"Measurement $\beta$ values")

for i in range(len(betas)):
    ax[1].bar(indices, generalized_Q[i],
              label = r"$beta = {:.2f}$".format(betas[i]))

ax[1].set_title("Population measurement statistics")
ax[1].set_xlabel("n")
ax[1].set_ylabel("Photon number probability")

hinton(rho, ax=ax[2])
ax[2].set_xlabel("Hinton plot of density matrix")
ax[1].legend()

plt.show()


